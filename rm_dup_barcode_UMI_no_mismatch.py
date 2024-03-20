
'''
This script accept a input sorted sam file, a output sam file, and a mismatch rate, then it will remove
duplicates based on the barcode + UMI (edit distance <= 1), and chromatin and start site, at the same
time, it will output the duplication number for each read, and generate the histogram plot for the read
per duplication number
'''
# using python3 updated 2024-03-12 used ChatGPT to convert so hopefully they did a good job

import sys
import numpy as np  
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
sys.setrecursionlimit(1000000)

def rm_dup_samfile(samfile, output_file, mismatch):
    with open(samfile) as f1, open(output_file, 'w') as f2, open(output_file+'.csv', 'w') as f3:
        pre_barcode = set()
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

                if start_site == pre_site and chrom_num == pre_chrom:
                    dup = False
                    if barcode_UMI not in pre_barcode:
                        pre_dup_num = cur_dup_num
                        cur_dup_num = 1
                        f2.write(line)
                        pre_barcode.add(barcode_UMI)
                        f3.write('%d\n' % pre_dup_num)
                    else:
                        cur_dup_num += 1
                else:
                    pre_dup_num = cur_dup_num
                    cur_dup_num = 1
                    f2.write(line)
                    pre_chrom = chrom_num
                    pre_site = start_site
                    pre_barcode = set()
                    pre_barcode.add(barcode_UMI)
                    if pre_dup_num != 0:
                        f3.write("%d\n" % pre_dup_num)

    
    '''
    #plot the histogram for the read duplication number
    dups = (pd.read_csv(output_file+'.csv', header=None))[0]
    fig = plt.figure()
    plt.hist(dups, bins=100)
    plt.xlabel("Duplication number")
    plt.ylabel("Read number")
    fig.savefig(output_file + '.png')
    '''


if __name__ == "__main__":
    samfile = sys.argv[1]
    output_file = sys.argv[2]
    mismatch = sys.argv[3]
    rm_dup_samfile(samfile, output_file, mismatch)
    
