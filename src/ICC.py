import networkx as nx
import pandas as pd
import numpy as np
import os
import sys, argparse

def get_parser():
    description = 'Acquisition of standardized networks'
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-d', '--ppi_dir', type = str, required = True, help = 'PPI network path')
    parser.add_argument('-s', '--species', type = str, required = True, help = 'Study species')
    parser.add_argument('-g', '--gene_group', type = str, required = True, help = 'Gene groups')
    parser.add_argument('-oh', '--one2one_Hs', type = str, required = True, help = '1:1 orthologos with human')
    parser.add_argument('-o', '--intermediate_o', type = str, required = True, help = 'ICCs for randomized network')
    return parser

def iccFun(A, B):
    """
    Interaction conservation coefficient (ICC)
    A: target set, for example human genes
    """
    intersec = list(set(A).intersection(set(B)))
    return len(intersec)/len(list(set(A)))

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

    for sp in species_dict.keys():
        if sp == args.species:
            stu_OG = {}
            with open(args.gene_group) as f:
                for row in f:
                    cols = row.strip().split('\t')
                    for gene in cols[1:]:
                        name = gene.split('|')
                        if name[0] == species_dict[sp]:
                            stu_OG[name[1]] = cols[0]
        else:   
            ref_ppi_file = args.ppi_dir + 'integrated_ppi_' + species_dict[sp] + '.txt'
            ref_G = nx.read_weighted_edgelist(ref_ppi_file, delimiter = '\t')

            ref_OG = {}
            with open(args.gene_group) as f:
                for row in f:
                    cols = row.strip().split('\t')
                    for gene in cols[1:]:
                        name = gene.split('|')
                        if name[0] == species_dict[sp]:
                            ref_OG[name[1]] = cols[0]
        
            for i in range(1000):
                rand_ppi_file = args.intermediate_o + species_dict[sp] + '/edge_list_' + str(i) + '.txt'
                rand_G = nx.read_weighted_edgelist(rand_ppi_file, delimiter = ' ')

                OG = []
                study = []
                species = []
                study_gene = []
                gene = []
                ICC = []
                with open(args.one2one_Hs) as f:
                    for row in f:
                        if row.startswith("OG"):
                            pass
                        else:
                            cols = row.strip().split('\t')
                            if cols[2] == sp:
                                if (cols[3] in list(rand_G.nodes())) & (cols[4] in list(ref_G.nodes())):
                                    a = [stu_OG[x] if x in stu_OG else x for x in list(rand_G.adj[cols[3]])]
                                    b = [ref_OG[x] if x in ref_OG else x for x in list(ref_G.adj[cols[4]])]
                                    OG.append(cols[0])
                                    study.append(cols[1])
                                    species.append(cols[2])
                                    study_gene.append(cols[3])
                                    gene.append(cols[4])
                                    ICC.append(iccFun(a, b))

                ortho_icc = pd.DataFrame({'OG': OG,
                                          'study': study,
                                          'species': species,
                                          'study_gene': study_gene,
                                          'gene': gene,
                                          'ICC': ICC})
                icc_out = args.intermediate_o + species_dict[sp] + '/ICC_edge_list_' + str(i) + '.txt'
                ortho_icc.to_csv(path_or_buf = icc_out, sep = '\t', index = False)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))









