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
    parser.add_argument('-o', '--intermediate_o', type = str, required = True, help = 'ICs for randomized network')
    return parser

def get_local_subgraph(G, node, radius = 2):
    """
    get a local subgraph for a node
    """
    nodes_within_radius = nx.single_source_shortest_path_length(G, node, cutoff = radius).keys()
    subgraph = G.subgraph(nodes_within_radius).copy()
    return subgraph

def isolated_components(G, node):
    """
    get isolated network components (ICs)
    """
    G_copy = G.copy()
    G_copy.remove_node(node)
    component_num = nx.number_connected_components(G_copy)
    k = G.degree(node)
    return round(component_num/k, 3)

def ic_output(ppi_f, ic_f):
    """
    ppi_f: input ppi network file
    in_f: output IC file
    """
    G = nx.read_weighted_edgelist(ppi_f, delimiter='\t')
    ICs = []
    for node in G.nodes():
        subG = get_local_subgraph(G, node, radius = 2)
        ic = isolated_components(subG, node)
        ICs.append(ic)
    IC_file = pd.DataFrame({'Node': G.nodes(),
                            'IC': ICs})
    IC_file.to_csv(path_or_buf = ic_f,
                   index = False, 
                   sep = '\t')

def ic_output2(ppi_f, ic_f):
    """
    ppi_f: input ppi network file
    in_f: output IC file
    """
    G = nx.read_weighted_edgelist(ppi_f, delimiter=' ')
    ICs = []
    for node in G.nodes():
        subG = get_local_subgraph(G, node, radius = 2)
        ic = isolated_components(subG, node)
        ICs.append(ic)
    IC_file = pd.DataFrame({'Node': G.nodes(),
                            'IC': ICs})
    IC_file.to_csv(path_or_buf = ic_f,
                   index = False, 
                   sep = '\t')

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
    # get ICs of nodes for each species' PPI network 
    ppi_files = [args.ppi_dir + 'integrated_ppi_' + sp + '.txt' for sp in species_dict.values() if sp != 'Hs']
    ic_files = [args.ppi_dir + 'orthologs_with_human/IC/IC_' + sp + '.txt' for sp in species_dict.values() if sp != 'Hs']
    combs = list(itertools.zip_longest(ppi_files, ic_files))
    with mp.Pool(6) as pool:
        pool.starmap(ic_output, combs)
        pool.close()
        pool.join()

    # get ICs of nodes for corresponding randomized human PPI network for each species
    for sp in species_dict.values():
        if sp != 'Hs':
            path =  args.intermediate_o + sp
            ppi_files = [path + "/edge_list_" + str(i) + ".txt" for i in range(1000)]
            ic_files = [path + "/IC_edge_list_" + str(i) + ".txt" for i in range(1000)]
            combs = list(itertools.zip_longest(ppi_files, ic_files))
            with mp.Pool(mp.cpu_count()-2) as pool:
                pool.starmap(ic_output2, combs)
                pool.close()
                pool.join()
    
if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))























