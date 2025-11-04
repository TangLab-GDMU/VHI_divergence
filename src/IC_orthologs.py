import networkx as nx
import pandas as pd
import numpy as np
import math
import itertools
import multiprocessing as mp
import random
import os
import sys, argparse

def get_parser():
    description = 'Acquisition of standardized networks'
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-oh', '--one2one_Hs', type = str, required = True, help = '1:1 orthologos with human')
    parser.add_argument('-o', '--intermediate_o', type = str, required = True, help = 'ICs for randomized network')
    parser.add_argument('-of', '--IC_out', type = str, required = True, help = 'out IC file')
    return parser

def hs_ICs(file, hs_gene):
    with open(file) as f:
        for row in f:
            if row.startswith("Node"):
                pass
            else:
                cols = row.strip().split('\t')
                if cols[0] == hs_gene:
                    return cols[1]

def calculate_z_score(ic, m, s):
    if s == 0:
        z = 0
    else:
        z = (ic - m) / s
    return z

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

    OG = []
    study = []
    species = []
    study_gene = []
    gene = []
    IC = []
    mean = []
    sd = []
    Z = []
    with open(args.one2one_Hs) as f:
        for row in f:
            if row.startswith("OG"):
                pass
            else:
                cols = row.strip().split('\t')
                sp = species_dict[cols[2]]

                # ICs of orthologs for human random network
                rand_ic_files = [args.intermediate_o + sp + "/IC_edge_list_" + str(i) + ".txt" for i in range(1000)]
                combs = list(itertools.product(rand_ic_files, [cols[3]]))
                with mp.Pool(mp.cpu_count()-10) as pool:
                    res = pool.starmap(hs_ICs, combs)
                    pool.close()
                    pool.join()
                res = [float(i) for i in res if i is not None]
                m, s = np.mean(res), np.std(res)

                # ICs of orthologs for orther species
                sp_IC_f = args.IC_out + 'IC_' + sp + '.txt'
                with open(sp_IC_f) as f2:
                    for row2 in f2:
                        if row2.startswith("Node"):
                            pass
                        else:
                            cols2 = row2.strip().split('\t')
                            if cols2[0] == cols[4]:
                                z = calculate_z_score(float(cols2[1]), m, s)
                                IC.append(cols2[1])
                                OG.append(cols[0])
                                study.append(cols[1])
                                species.append(cols[2])
                                study_gene.append(cols[3])
                                gene.append(cols[4])
                                mean.append(m)
                                sd.append(s)
                                Z.append(z)

    ortho_ic = pd.DataFrame({'OG': OG,
                             'study': study,
                             'species': species,
                             'study_gene': study_gene,
                             'gene': gene,
                             'IC': IC,
                             'mean': mean,
                             'sd': sd,
                             'Z': Z})
    ortho_ic_out = args.IC_out + 'IC_ortho1to1_with_human.txt'
    ortho_ic.to_csv(path_or_buf = ortho_ic_out, sep = '\t', index = False)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))

