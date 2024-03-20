
'''
This script accept a input sorted sam file, a output sam file, and a mismatch rate, then it will remove
duplicates based on the barcode + UMI (edit distance <= 1), and chromatin and start site, at the same
time, it will output the duplication number for each read, and generate the histogram plot for the read
per duplication number
'''
from Levenshtein import distance
import sys
import matplotlib as mpl
import numpy as np
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

sys.setrecursionlimit(1000000)

def rm_dup_samfile(samfile, output_file, mismatch):
    with open(samfile, 'r') as f1, open(output_file, 'w') as f2, open(output_file+'.csv', 'w') as f3:
        pre_barcode = []
        pre_line = []
        unique_id = []
        pre_chrom = 0
        pre_site = 0
        dup = False
        mismatch = int(mismatch)

        pre_dup_num = 0
        cur_dup_num = 0

        for line in f1:
            if line[0] == '@':
                f2.write(line)
            else:
                name = (((line.split('\t'))[0]).split(','))
                barcode_UMI = name[0] + name[1]
                chrom_num = (line.split('\t'))[2]
                start_site = (line.split('\t'))[3]

                if (start_site == pre_site) and (chrom_num == pre_chrom):
                    if barcode_UMI in pre_barcode:
                        cur_dup_num += 1
                    else:
                        pre_barcode.append(barcode_UMI)
                        pre_line.append(line)

                else:
                    if pre_barcode:
                        unique_id = index_unique(pre_barcode, mismatch)
                        for i in unique_id:
                            f2.write(pre_line[i])
                        cur_dup_num = cur_dup_num + len(pre_barcode) - len(unique_id)

                    pre_dup_num = cur_dup_num
                    cur_dup_num = 1
                    unique_id = []
                    pre_barcode = []
                    pre_line = []
                    pre_chrom = chrom_num
                    pre_site = start_site

                    pre_barcode.append(barcode_UMI)
                    pre_line.append(line)
                    if pre_dup_num != 0:
                        f3.write(str(pre_dup_num))
                        f3.write('\n')

        # also count for the last reads
        if pre_barcode:
            unique_id = index_unique(pre_barcode, mismatch)
            for i in unique_id:
                f2.write(pre_line[i])
            cur_dup_num = cur_dup_num + len(pre_barcode) - len(unique_id)

def index_unique(barcodes_list, cutoff_value):
    list_length = len(barcodes_list)
    distance_matrix = np.arange(list_length * list_length).reshape(list_length, list_length)

    for i in range(list_length):
        for j in range(list_length):
            if j < i:
                distance_matrix[i][j] = distance_matrix[j][i]
            elif j == i:
                distance_matrix[i][j] = 0
            else:
                distance_matrix[i][j] = distance(barcodes_list[i], barcodes_list[j])

    unique_index = []
    non_visited_indexes = list(range(list_length))
    for i in range(list_length):
        if i in non_visited_indexes:
            unique_index.append(i)
            non_visited_indexes = rm_dup(non_visited_indexes, distance_matrix, i, cutoff_value)

    return unique_index

def rm_dup(indexes, distance_matrix, n, cutoff_value):
    if n in indexes:
        indexes.remove(n)

    all_indexes = indexes[:]
    for i in range(len(all_indexes)):
        if all_indexes[i] in indexes:
            if distance_matrix[all_indexes[i], n] <= cutoff_value:
                indexes = rm_dup(indexes, distance_matrix, all_indexes[i], cutoff_value)
    return indexes

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 script.py <samfile> <output_file> <mismatch>")
        sys.exit(1)

    samfile = sys.argv[1]
    output_file = sys.argv[2]
    mismatch = sys.argv[3]
    rm_dup_samfile(samfile, output_file, mismatch)
    
