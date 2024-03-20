
# This script take as input a gtf file, a sam file input folder, a sample ID file, and the number of cores, then run
# the count exons and genes in parallel, and then output the report file and output file into the output file

#converted to python3 by ChatGPT

import itertools
import collections
import numpy as np
import pandas as pd
from multiprocessing import Pool
from multiprocessing import *
import HTSeq
import sys
from functools import partial
import logging

def sciRNAseq_count(sample, input_folder, exons, genes, gene_end, gene_annotat, sample_ID):
    input_sam = f"{input_folder}/{sample}.sam"
    report = f"{input_folder}/{sample}.report"
    count_output = f"{input_folder}/{sample}.count"

    counts = collections.Counter()
    sam_file = input_sam
    almnt_file = HTSeq.SAM_Reader(sam_file)
    sam_name = sample
    cell_ID = sample_ID.index(sample) + 1

    perfect_inter_exon = 0
    nearest_inter_exon = 0
    perfect_combine_exon = 0
    nearest_combine_exon = 0
    perfect_inter_gene = 0
    nearest_inter_gene = 0
    perfect_combine_gene = 0
    nearest_combine_gene = 0
    
    print(f"Start read the input file: {sam_file}....")

    for alnmt in almnt_file:
        if not alnmt.aligned:
            counts["_unmapped"] += 1
            continue
        
        if alnmt.iv.chrom not in genes.chrom_vectors:
            counts["_unmapped"] += 1
            continue

        gene_id_intersect = set()
        gene_id_combine = set()
        inter_count = 0
        for cigop in alnmt.cigar:
            if cigop.type != "M":
                continue

            for iv,val in exons[cigop.ref_iv].steps():
                gene_id_combine |= val
                if inter_count == 0:
                    gene_id_intersect |= val
                    inter_count += 1
                else:
                    gene_id_intersect &= val
        
        if len(gene_id_intersect) == 1:
            gene_id = list(gene_id_intersect)[0]
            counts[gene_id] += 1
            perfect_inter_exon += 1
        elif len(gene_id_intersect) > 1:
            gene_id = find_nearest_gene(alnmt.iv.end_d, gene_id_intersect, gene_end)
            counts[gene_id] += 1
            nearest_inter_exon += 1
        else:
            if len(gene_id_combine) == 1:
                gene_id = list(gene_id_combine)[0]
                counts[gene_id] += 1
                perfect_combine_exon += 1
            elif len(gene_id_combine) > 1:
                gene_id = find_nearest_gene(alnmt.iv.end_d, gene_id_combine, gene_end)
                counts[gene_id] += 1
                nearest_combine_exon += 1
            else:
                gene_id_intersect = set()
                gene_id_combine = set()
                inter_count = 0
                for cigop in alnmt.cigar:
                    if cigop.type != "M":
                        continue
                    for iv,val in genes[cigop.ref_iv].steps():
                        gene_id_combine |= val
                        if inter_count == 0:
                            gene_id_intersect |= val
                            inter_count += 1
                        else:
                            gene_id_intersect &= val

                if len(gene_id_intersect) == 1:
                    gene_id = list(gene_id_intersect)[0] + "_intron"
                    counts[gene_id] += 1
                    perfect_inter_gene += 1

                elif len(gene_id_intersect) > 1:
                    gene_id = find_nearest_gene(alnmt.iv.end_d, gene_id_intersect, gene_end) + "_intron"
                    counts[gene_id] += 1
                    nearest_inter_gene += 1

                else:
                    if len(gene_id_combine) == 1:
                        gene_id = list(gene_id_combine)[0] + "_intron"
                        counts[gene_id] += 1
                        perfect_combine_gene += 1

                    elif len(gene_id_combine) > 1:
                        gene_id = find_nearest_gene(alnmt.iv.end_d, gene_id_combine, gene_end) + "_intron"
                        counts[gene_id] += 1
                        nearest_combine_gene += 1

                    else:
                        counts["_no_feature"] += 1
                        
    print("File name: ", sam_file)
    print("1: Perfect intersect exon match: ", perfect_inter_exon)
    print("2: Nearest intersect exon match: ", nearest_inter_exon)
    print("3: Perfect combine exon match: ", perfect_combine_exon)
    print("4: Nearest combine exon match: ", nearest_combine_exon)
    print("5: Perfect intersect gene match: ", perfect_inter_gene)
    print("6: Nearest intersect gene match: ", nearest_inter_gene)
    print("7: Perfect combine gene match: ", perfect_combine_gene)
    print("8: Nearest combine gene match: ", nearest_combine_gene)
    print("9: ambiguous match for exons: ", counts["_ambiguous"])
    print("10: ambiguous match for genes: ", counts["_ambiguous_intron"])
    print("11: No match: ",  counts["_no_feature"])
    print("Sam file analysis finished~")
    
    with open(report, 'w') as report:
        report.write(f"1,{cell_ID},{perfect_inter_exon}\n")
        report.write(f"2,{cell_ID},{nearest_inter_exon}\n")
        report.write(f"3,{cell_ID},{perfect_combine_exon}\n")
        report.write(f"4,{cell_ID},{nearest_combine_exon}\n")
        report.write(f"5,{cell_ID},{perfect_inter_gene}\n")
        report.write(f"6,{cell_ID},{nearest_inter_gene}\n")
        report.write(f"7,{cell_ID},{perfect_combine_gene}\n")
        report.write(f"8,{cell_ID},{nearest_combine_gene}\n")
        report.write(f"9,{cell_ID},{counts['_ambiguous']}\n")
        report.write(f"10,{cell_ID},{counts['_ambiguous_intron']}\n")
        report.write(f"11,{cell_ID},{counts['_no_feature']}\n")
        
    with open(count_output, 'w') as count_output:
        for gene in counts:
            if (gene in ["_unmapped", "_ambiguous", "_ambiguous_intron", "_no_feature"]):
                continue
            else:
                line = f"{gene_annotat.loc[gene,4]},{cell_ID},{counts[gene]}\n"
                count_output.write(line)
    return 0

def find_nearest_gene(al_end, gene_id_intersect, gene_end):
    gene_id_end = {}
    for gene in gene_id_intersect:
        if gene in gene_end:
            gene_id_end[gene] = (abs(np.array(list(gene_end[gene])) - al_end)).min()
        else:
            print("****************Found one gene without transcript annotation*****************", "Gene name: ", gene)
    
    gene_end_min = np.min(list(gene_id_end.values()))
    count = 0
    for gene in gene_id_end:
        if (gene_id_end[gene] < gene_end_min + 100):
            count += 1
            gene_id = gene
    if count > 1:
        gene_id = "_ambiguous"
    
    return gene_id

def sciRNA_count_parallel(gtf_file, input_folder, sample_ID, core_number):
    gtf_file = HTSeq.GFF_Reader(gtf_file, end_included=True)
    gene_annotat_file = f"{input_folder}/gene_name_annotate.txt"
    cell_annotat_file = f"{input_folder}/cell_annotate.txt"
    report_annotate_file = f"{input_folder}/report_annotate.txt"
    
    gene_annotat = open(gene_annotat_file, "w")
    cell_annotat = open(cell_annotat_file, "w")
    report_annotate = open(report_annotate_file, "w")
    
    exons = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    genes = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    gene_end = {}
    exon_n = 0
    gene_n = 0
    transcript_n = 0
    gene_count = 0
    print("Start generating exon genomic arrays....")
    print("Start generating gene genomic arrays....")
    print("Start calculating transcript end of genes....")

    for feature in gtf_file:
        if feature.type == "exon":
            exon_n += 1
            exons[ feature.iv ] += feature.attr["gene_id"]
        elif feature.type == "gene":
            gene_n +=1
            genes[ feature.iv ] += feature.attr["gene_id"]
            gene_count += 1
            message = f"{feature.attr['gene_id']},{feature.attr['gene_type']},exon,{feature.attr['gene_name']},{gene_count}\n"
            gene_annotat.write(message)
            gene_count += 1
            message = f"{feature.attr['gene_id']}_intron,{feature.attr['gene_type']},intron,{feature.attr['gene_name']}_intron,{gene_count}\n"
            gene_annotat.write(message)
        elif feature.type == "transcript":
            transcript_n += 1
            if feature.attr["gene_id"] in gene_end.keys():
                gene_end[ feature.attr["gene_id"] ].add(feature.iv.end_d)
            else:
                gene_end[ feature.attr["gene_id"] ] = set()
                gene_end[ feature.attr["gene_id"] ].add(feature.iv.end_d)

    print("Detected gene number: ", gene_n)
    print("Detected transcript number: ", transcript_n)
    print("Detected exon number: ", exon_n)
    
    gene_annotat.close()
    
    gene_annotat = pd.read_csv(gene_annotat_file, header=None)
    gene_annotat.index =  gene_annotat[0]
    
    sample_ID = list(pd.read_csv(sample_ID, header=None)[0])
    
    cell_count = 0
    for i in sample_ID:
        cell_count += 1
        message = f"{i},{cell_count}\n"
        cell_annotat.write(message)
    cell_annotat.close()
    
    report_annotate.write("1, Perfect intersect exon match\n")
    report_annotate.write("2, Nearest intersect exon match\n")
    report_annotate.write("3, Perfect combine exon match\n")
    report_annotate.write("4, Nearest combine exon match\n")
    report_annotate.write("5, Perfect intersect gene match\n")
    report_annotate.write("6, Nearest intersect gene match\n")
    report_annotate.write("7, Perfect combine gene match\n")
    report_annotate.write("8, Nearest combine gene match\n")
    report_annotate.write("9, ambiguous match for exons\n")
    report_annotate.write("10, ambiguous match for genes\n")
    report_annotate.write("11, No match\n")
    report_annotate.close()
    
    p = Pool(processes = int(core_number))
    func = partial(sciRNAseq_count, input_folder=input_folder, exons=exons, genes=genes, gene_end=gene_end, gene_annotat=gene_annotat, sample_ID=sample_ID)
    result = p.map(func, sample_ID)
    p.close()
    p.join()

    print("All analysis done~")

if __name__ == "__main__":
    gtf_file = sys.argv[1]
    input_folder = sys.argv[2]
    sample_ID = sys.argv[3]
    core_number = sys.argv[4]
    sciRNA_count_parallel(gtf_file, input_folder, sample_ID, core_number)
    
