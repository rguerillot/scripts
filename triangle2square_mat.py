#!/usr/bin/python3

import sys
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Convert triangle phylip distance matrix to square (full phylip distance matrix')
parser.add_argument('triangle', type=str)
parser.add_argument('square', type=str)
parser.add_argument("-c", "--csv", help = "output a csv table wich column header instead of phylip", required = False, default = False)

args = parser.parse_args()
print("in file: " + args.triangle)
print("out file: " + args.square)

distances = []

def lower_triangle_to_full_matrix(filename, fileout):
    num_lines_in_file = sum(1 for line in open(filename))
    distances = []
    sample_names = []
    sample_nb=0

    with open(filename) as f:
        next(f) # skip sample count line
        for line in f:
            sample_nb+=1
            print('converting line: ' + str(sample_nb) + "/" + str(num_lines_in_file-1) , end='\r')
            elements = line.strip().split('\t')
            sample_names.append(elements[0])
            row = [float(e) for e in elements[1:]]
            row.extend([0.0] * (num_lines_in_file-1-len(row)))
            distances.append(row)
        np_array = np.asarray(distances)
        index_upper = np.triu_indices(num_lines_in_file-1)
        np_array[index_upper] = np_array.T[index_upper]
    if args.csv:
        return pd.DataFrame(np_array, columns=sample_names, index=sample_names).to_csv(fileout)
    else:
        text_file = open(fileout, "w")
        text_file.write(str(sample_nb) + "\n")
        text_file.close()
        return pd.DataFrame(np_array, columns=sample_names, index=sample_names).to_csv(fileout, sep='\t', header = False, mode='a')

lower_triangle_to_full_matrix(args.triangle, args.square)
