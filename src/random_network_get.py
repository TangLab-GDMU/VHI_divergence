import networkx as nx
import pandas as pd
import numpy as np
import math
import itertools
import multiprocessing as mp
import random
import os
import sys, argparse
from scipy.optimize import curve_fit

def get_parser():
    description = 'Acquisition of standardized networks'
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-d', '--ppi_dir', type = str, required = True, help = 'PPI network path')
    parser.add_argument('-s', '--species', type = str, required = True, help = 'Study species')
    parser.add_argument('-se', '--seed_set', type = str, required = True, help = 'Seed set')
    parser.add_argument('-o', '--output', type = str, required = True, help = 'Out put path')
    return parser

def net_size_estimation(stu_N, ref, ppi_dir):
    """
    stu_N: number of edge selected for study species
    ref: network file for reference species
    ppi_dir: network file path
    IRR = -0.215
    output: node number and edge number
    """
    ppi_f = ppi_dir + 'integrated_ppi_' + ref + '.txt'
    ref_G = nx.read_weighted_edgelist(ppi_f, delimiter='\t')
    N = len(list(ref_G.nodes))
    M = len(list(ref_G.edges))
    p = M/math.comb(N, 2)
    A = N * p
    if ref == 'Xl':    
        stu_A = (2**0.215) * A
        stu_p = stu_A/stu_N
        stu_M = math.comb(stu_N, 2) * stu_p
        return stu_N, stu_M
    else:
        stu_A = (2**-0.215) * A
        stu_p = stu_A/stu_N
        stu_M = math.comb(stu_N, 2) * stu_p
        return stu_N, stu_M

def weighted_edge_sampler(G, num_edges, seed=2914):
    """
    weight for each edge (e = (u,v)): w = degree(u) + degree(v)
    """
    # calculating the weights for each edge
    edge_weights = []
    for u,v in G.edges():
        w = G.degree(u) + G.degree(v)
        edge_weights.append(w)
    edges = list(G.edges())
    # converting to probability distribution
    probs = np.array(edge_weights, dtype = float)
    probs /= probs.sum()
    # sampling with replacement
    np.random.seed(seed)
    chosen_idx = np.random.choice(len(edges), size = num_edges, replace = False, p = probs)
    return [edges[i] for i in chosen_idx]

def remove_low_degree_nodes(G, min_degree = 1):
    nodes_to_remove = [n for n,d in G.degree() if d < min_degree]
    G_pruned = G.copy()
    G_pruned.remove_nodes_from(nodes_to_remove)
    return G_pruned

def sample_subgraph(G, num_edges, seed=2914):
    nodes = list(G.nodes())
    # 1) sampling
    sampled_edges = weighted_edge_sampler(G, num_edges, seed)
    # 2) constructing subgraph
    H = nx.Graph()
    H.add_nodes_from(nodes)
    H.add_edges_from(sampled_edges)
    H2 = remove_low_degree_nodes(H)
    return H2

def get_node_edge_num(stu, ppi_dir, species_dict, seed=2914):
    stu_f = ppi_dir + 'integrated_ppi_' + species_dict[stu] + '.txt'
    stu_G = nx.read_weighted_edgelist(stu_f, delimiter='\t')
    sp_size = {}
    for sp in species_dict.values():
        if sp != 'Hs':
            est_edge_num = net_size_estimation(len(stu_G.nodes()),sp, ppi_dir)
            for i in range(int(est_edge_num[1]), 5000, -100):
                stu_H = sample_subgraph(stu_G, i, seed)
                est_edge_num2 = net_size_estimation(len(stu_H.nodes()), sp, ppi_dir)
                if len(stu_H.edges()) < est_edge_num2[1]:
                    #print(sp, "\t", len(stu_H.edges()), "\t", est_edge_num2[0], "\t", est_edge_num2[1])
                    sp_size[sp] = [len(stu_H.edges()), est_edge_num2[0], est_edge_num2[1]]
                    break
    return sp_size

def run(args):
    species_dict = {
        'H. sapiens': 'Hs',
        'M. musculus': 'Mm',
        'R. norvegicus': 'Rn',
        'X. laevis': 'Xl',
        'D. melanogaster': 'Dm',
        'C. elegans': 'Ce',
        'S. cerevisiae': 'Sc'
    }

    stu_f = args.ppi_dir + 'integrated_ppi_' + species_dict[args.species] + '.txt'
    stu_G = nx.read_weighted_edgelist(stu_f, delimiter='\t')
    species_size = get_node_edge_num(args.species, args.ppi_dir, species_dict)

    for sp in species_dict.values():
        if sp != 'Hs':
            os.makedirs(args.output + sp, exist_ok = True)
            edge_size = int(species_size[sp][2])
    
            for i in range(1000):
                np.random.seed(i)
                H2 = sample_subgraph(stu_G, edge_size, i)
                print(sp, "\t", len(H2.nodes()), "\t", len(H2.edges()))
                f = args.output + sp + "/edge_list_" + str(i) + ".txt"
                nx.write_edgelist(H2, f, data = False)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))












