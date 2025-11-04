import os
import sys, argparse
from Bio import SeqIO

def get_parser():
    description = 'Acquisition of sequences'
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-oh', '--one2one_Hs', type = str, required = True, help = '1:1 orthologos with human')
    parser.add_argument('-ifa', '--input_fa', type = str, required = True, help = 'fasta for all genes')
    parser.add_argument('-ofa', '--output_fa', type = str, required = True, help = 'fasta for pairwise genes')
    return parser

def extract_sequence(fasta_file, target_id, output_file=None):
    """sequence extraction of specified IDs"""
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == target_id:
            if output_file:
                with open(output_file, "a") as handle:
                    SeqIO.write(record, handle, "fasta")
            return record
    return None

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

    with open(args.one2one_Hs) as f:
        for row in f:
            if row.startswith("OG"):
                pass
            else:
                cols = row.strip().split('\t')
                stu_sp = species_dict[cols[1]]
                sp = species_dict[cols[2]]
                stu_seq_f = args.input_fa + stu_sp + '.cds.fa'
                sp_seq_f = args.input_fa + sp + '.cds.fa'
                stu_seq = extract_sequence(stu_seq_f, cols[3])
                sp_seq = extract_sequence(sp_seq_f, cols[4])
                #print(cols[3])
                #print(cols[4])
                if stu_seq and sp_seq:
                    out_f = args.output_fa + stu_sp + "_" + sp + "_" + cols[3] + "_" + cols[4] + ".fasta"
                    #print(out_f)
                    extract_sequence(stu_seq_f, cols[3], out_f)
                    extract_sequence(sp_seq_f, cols[4], out_f)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))

