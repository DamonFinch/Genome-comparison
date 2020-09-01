from CHROMEISTER import *
from common_functions import *
import sys
import os
import numpy as np
import pandas as pd
import csv
from main import *
import argparse

import resource

z1 = [0, 4, 8, 12, 16, 20, 24, 28]
z2 = [0, 1, 8, 12, 14, 15, 16, 17, 24, 28, 30, 31]
z3 = [3, 5, 7, 11, 13, 15, 19, 21, 23, 27, 29, 31]
z4 = [2, 5, 6, 10, 13, 14, 18, 21, 22, 26, 29, 30]
z5 = [12, 16, 20, 28]
z6 = [12, 28]
z7 = [28]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='CHROMEISTER')
    parser.add_argument('--path', required=True, type=str, help='path to fastas')
    parser.add_argument('--key-len', type=int, default=12, help='length of key for hit')
    parser.add_argument('--k-mer-len', type=int, default=32, help='length of k-mers')
    parser.add_argument('--dim', type=int, default=1000, help='dimensions of dotplots')
    parser.add_argument("--verbose", help="increase output verbosity", action="store_true")
    args = parser.parse_args()

    path = args.path
    key_len = args.key_len
    kmer_len = args.k_mer_len
    dim = args.dim
    log = args.verbose

    all_files = [f for f in os.listdir(path) if ((not f.startswith('.')) and (os.path.isfile(os.path.join(path, f))) )]    
    mat = np.zeros((len(all_files), len(all_files)))
    dic_scores = {}
    for i, file1 in enumerate(all_files):
        for j in range(i + 1 , len(all_files)):
            file2 = all_files[j]
            score, query_name, ref_name = chromeister(path , file1, file2, z1, dim=dim, log=log)
            mat[j][i] = score
            if query_name not in dic_scores:
                dic_scores[query_name] = {}
            dic_scores[query_name][ref_name] = mat[j][i]
            print("---------------------------------------------------------------------")
    df = pd.DataFrame(dic_scores)
    df.to_csv(path + '/results/' + 'distmat.csv')


    # print(mat)

print("Memory Usage :" + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000) + str(' KB'))