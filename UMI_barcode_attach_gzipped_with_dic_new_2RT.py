"""
Created on Tue Apr  5 22:16:54 2016

@author: Junyue
beth tweaking this to add pcr well to header also, this one is not quite done yet, I haven't worked on it in a while

"""

import subprocess
import sys
from Levenshtein import distance
import gzip
from multiprocessing import Pool
from functools import partial
import pickle
import csv #new
import os #new

'''
    this script accept a read1 file, a read2 file, a output_file, a ligation barcode list,
    a oligodT barcode list,
    and mismatch rate, then it open the read1 and read2, output file,
    then extract the barcode and UMI sequence in the read 1 file, and convert the
    barcode to the real barcode in the barcode list based on the mismatch rate,
    then it attach the barcode and UMI sequence to the read name of the read2 file
'''    
    
def UMI_attach_read2_barcode_list(pcrwell, input_folder, output_folder, ligation_barcode_list, RT_barcode_list, RTsamplename_dict, RTprimer_dict, RTwells_dict, LIGwells_dict, mismatch_rate = 1):
    #open the read1, read2, and output file
    Read1 = input_folder + "/" + pcrwell + ".R1.fastq.gz"
    Read2 = input_folder + "/" + pcrwell + ".R2.fastq.gz"
    output_file = output_folder + "/" + pcrwell + ".R2.fastq.gz"
    mismatch_rate = int(mismatch_rate)
    f1 = gzip.open(Read1)
    f2 = gzip.open(Read2)
    f3 = gzip.open(output_file, 'wb')
    #suffix = ".fastq.gz"
        
    line1 = f1.readline()
    line2 = f2.readline()
    total_line = 0
    filtered_line = 0
    
    while (line1):
        total_line += 1
        line1 = f1.readline()
        #print("read1: ", line1)
        # first check if the ligation barcode match with the barcode
        tmp_lig = line1[0:10]

        #print("ligation barcode: ", ligation_bc_match)
        if tmp_lig in ligation_barcode_list:
            ligation_bc_match = ligation_barcode_list[tmp_lig]
            if ligation_bc_match in LIGwells_dict:
                ligindx = LIGwells_dict[ligation_bc_match]
                # check RT barcode
                target_RT = line1[len(ligation_bc_match) + 14 : len(ligation_bc_match) + 24]
                #print("target_RT: ", target_RT)

                if target_RT in RT_barcode_list:
                    barcode = RT_barcode_list[target_RT]
                    if barcode in RTsamplename_dict: #now check if it's even an RT barcode that you used and connect to samplename and RTprimer name
                        samplename = RTsamplename_dict[barcode]
                        RTprimer = RTprimer_dict[barcode]
                        RTwell = RTwells_dict[barcode]
                        filtered_line += 1
                        UMI = line1[len(ligation_bc_match) + 6 : len(ligation_bc_match) + 14]
                        first_line = '@LIG-' + ligindx + '_RT-' + RTwell + ',' + UMI + ',' + samplename + ',' + RTprimer + '\n'  #adding samplename and RTprimername
                        #print("read2: ", first_line)
                        f3.write(first_line)

                        second_line = f2.readline()
                        f3.write(second_line)

                        third_line = f2.readline()
                        f3.write(third_line)

                        four_line = f2.readline()
                        f3.write(four_line)

                        line2 = f2.readline()
            
                        # output_file = output_folder + "/" + samplename + ".fastq.gz"
                        # with gzip.open(output_file, 'ab') as f3: 
                            # filtered_line += 1
                            # UMI = line1[len(ligation_bc_match) + 6 : len(ligation_bc_match) + 14]
                            # first_line = '@LIG-' + ligation_bc_match + '_RT-' + barcode + '_PCR-' + pcrwell + ',' + UMI + ',' + samplename + ',' + RTprimer + '\n'  #adding pcr well (pcrwell) here and samplename
                            # #print("read2: ", first_line)
                            # f3.write(first_line)

                            # second_line = f2.readline()
                            # f3.write(second_line)

                            # third_line = f2.readline()
                            # f3.write(third_line)
                            # four_line = f2.readline()
                            # f3.write(four_line)

                            # line2 = f2.readline()

                    else:
                        line2 = f2.readline()
                        line2 = f2.readline()
                        line2 = f2.readline()
                        line2 = f2.readline()                 
                else:
                    line2 = f2.readline()
                    line2 = f2.readline()
                    line2 = f2.readline()
                    line2 = f2.readline()
                    
            else:
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                
        else:
            line2 = f2.readline()
            line2 = f2.readline()
            line2 = f2.readline()
            line2 = f2.readline()

        line1 = f1.readline()
        line1 = f1.readline()
        line1 = f1.readline()

    f1.close()
    f2.close()
    #f3.close()
    
    
    print("PCR well: %s, total line: %f, filtered line: %f, filter rate: %f" 
          %(pcrwell, total_line, filtered_line, float(filtered_line) / float(total_line)))

# this function accept an input folder and a output folder and then generate the output file with the index
def attach_UMI_files(input_folder, pcrwell, output_folder, ligation_barcode_file, LIGlist, RT_barcode_file, RT_sample, core):
    
    init_message = '''
    --------------------------start attaching UMI-----------------------------
    input folder: %s
    PCR well list: %s
    output_folder: %s
    ligation barcode file: %s
    ligation list: %s
    RT barcode file: %s
    RT samplesheet: %s
    ______________________________________________________________________________________________________________________________________________________
    ''' %(input_folder, pcrwell, output_folder, ligation_barcode_file, LIGlist, RT_barcode_file, RT_sample)
    
    print(init_message)
    
    print("Load ligation barcode dictionary...")
    
    # generate the ligation barcode list
    barcodes = open(ligation_barcode_file, "rb")
    ligation_barcode_list = pickle.load(barcodes)
    barcodes.close()
    
    print("Load RT barcode dictionary...")
    
    # generate the RT barcode list:
    barcodes = open(RT_barcode_file, "rb")
    RT_barcode_list = pickle.load(barcodes)
    #print(RT_barcode_list)
    barcodes.close()
    
    #new
    print("Loading sample names...")
    RTsamplenames = open(RT_sample, "r")
    RTreader = csv.DictReader(RTsamplenames)
    RTsamplename_dict = {rows["RTindex"]:rows["SampleName"] for rows in RTreader}
    samplelist = list(set(val for dic in RTsamplename_dict for val in RTsamplename_dict.values()))
    #print(RTsamplename_dict)
    print(samplelist)
    print("Loading RT primers")
    #adding this for 2RT primer info
    RTprimernames = open(RT_sample, "r")
    RTreader2 = csv.DictReader(RTprimernames)
    RTprimer_dict = {rows["RTindex"]:rows["RTprimer"] for rows in RTreader2}
    RTprimerlist = list(set(val for dic in RTprimer_dict for val in RTprimer_dict.values()))
    #print(RTprimer_dict)

    RTwells = open(RT_sample, "r")
    RTreader3 = csv.DictReader(RTwells)
    RTwells_dict = {rows["RTindex"]:rows["RTwell"] for rows in RTreader3}
    RTwellslist = list(set(val for dic in RTwells_dict for val in RTwells_dict.values()))
    #print(RTwells_dict)
    print(RTprimerlist)
    
    LIGwells = open(LIGlist, "r")
    RTreader4 = csv.DictReader(LIGwells)
    LIGwells_dict = {rows["LIGindex"]:rows["LIGwell"] for rows in RTreader4}
    LIGwellslist = list(set(val for dic in LIGwells_dict for val in LIGwells_dict.values()))
    #print(LIGwellslist)
    
    #for each item in the pcrwell list, use the read1 file, read2 file, output file
    # and barcode_list to run UMI_attach_read2_barcode_list
    pcrwell_file = open(pcrwell)
    pcrwell_list = []
    for line in pcrwell_file:
        pcrwell = line.strip()
        pcrwell_list.append(pcrwell)
    pcrwell_file.close()
    
    # parallele for the functions
    p = Pool(processes = int(core))
    #print("Processing core number: ", core_number)
    func = partial(UMI_attach_read2_barcode_list, input_folder = input_folder, output_folder=output_folder, ligation_barcode_list = ligation_barcode_list, RT_barcode_list=RT_barcode_list, RTsamplename_dict = RTsamplename_dict, RTprimer_dict = RTprimer_dict, RTwells_dict = RTwells_dict, LIGwells_dict= LIGwells_dict, mismatch_rate = 1)
    #sciRNAseq_count(pcrwell, input_folder, exons, genes, gene_end)
    result = p.map(func, pcrwell_list)
    p.close()
    p.join()
    
    #print the completion message
    com_message = '''~~~~~~~~~~~~~~~UMI attachment done~~~~~~~~~~~~~~~~~~'''
    print(com_message)
    
if __name__ == "__main__":
    input_folder = sys.argv[1]
    pcrwell = sys.argv[2]
    output_folder = sys.argv[3]
    ligation_barcode_file = sys.argv[4]
    LIGlist = sys.argv[5] #you have to get rid of this, I don't think you'll use it
    RT_barcode_file = sys.argv[6]
    RT_sample = sys.argv[7]
    core=sys.argv[8]
    attach_UMI_files(input_folder, pcrwell, output_folder, ligation_barcode_file, LIGlist, RT_barcode_file, RT_sample, core)
    
