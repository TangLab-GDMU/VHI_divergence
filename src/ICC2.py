import pandas as pd
import os
import sys, argparse

def get_parser():
    description = 'Acquisition of standardized networks'
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-i', '--intermediate_d', type = str, required = True, help = 'ICCs for randomized network')
    parser.add_argument('-s', '--species', type = str, required = True, help = 'Study species')
    parser.add_argument('-o', '--ICC_stat', type = str, required = True, help = 'ICC results')
    return parser

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
        if sp != args.species:
            icc_df = pd.DataFrame(columns = ['OG', 'study', 'species', 'study_gene', 'gene', 'ICC'])
            for i in range(1000):
                file = args.intermediate_d + species_dict[sp] + '/ICC_edge_list_' + str(i) + '.txt'
                icc = pd.read_csv(file, sep="\t")
                icc_df = pd.concat([icc_df, icc])

            result = icc_df.groupby(['OG', 'study', 'species', 'study_gene', 'gene']).agg({
                'ICC':['mean', 'std'],
                'ICC':['mean', 'std']
            })
            result.reset_index(inplace = True)
            result.columns = ['OG', 'study', 'species', 'study_gene', 'gene', 'mean', 'std']
            out_f = args.ICC_stat + 'ICC_' + species_dict[sp] + '.txt'
            result.to_csv(path_or_buf = out_f,
                          index = False,
                          sep = '\t')

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))