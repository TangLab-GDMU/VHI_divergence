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
    description = 'Acquisition of degree exponent'
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-f', '--ppi_file', type = str, required = True, help = 'PPI network')
    parser.add_argument('-d', '--delimiter', type = str, required = False, help = 'delimiter')
    return parser

def power_law(x, gamma, c):
    return c * (x ** (-gamma))

def get_net_from_file(f, delimiter = '\t'):
    G = nx.read_weighted_edgelist(f, delimiter='\t')
    return G

def get_degree_exponent(G):
    degrees = [d for _, d in G.degree()]
    k_min = min(degrees)
    k_max = max(degrees)
    hist, bin_edges = np.histogram(degrees, 
                               bins = np.logspace(np.log10(k_min), np.log10(k_max),50),
                              density = True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    popt, _ = curve_fit(power_law, bin_centers[hist > 0],
                        hist[hist > 0],
                        p0 = [2.0, 1.0])
    gamma_estimated = popt[0]
    return round(gamma_estimated, 3)

def run(args):
    G = get_net_from_file(args.ppi_file)
    gamma = get_degree_exponent(G)
    print(gamma)


if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))