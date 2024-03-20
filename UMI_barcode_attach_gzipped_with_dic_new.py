#chatgpt conversion to python3

import sys
import gzip
import pickle
import csv
from multiprocessing import Pool
from functools import partial

def UMI_attach_read2_barcode_list(pcrwell, input_folder, output_folder, ligation_barcode_list, RT_barcode_list, RTsamplename_dict, RTwells_dict, LIGwells_dict, mismatch_rate=1):
    Read1 = input_folder + "/" + pcrwell + ".R1.fastq.gz"
    Read2 = input_folder + "/" + pcrwell + ".R2.fastq.gz"
    output_file = output_folder + "/" + pcrwell + ".R2.fastq.gz"
    mismatch_rate = int(mismatch_rate)

    with gzip.open(Read1, 'rt') as f1, gzip.open(Read2, 'rt') as f2, gzip.open(output_file, 'wt') as f3:
        line1 = f1.readline()
        line2 = f2.readline()
        total_line = 0
        filtered_line = 0

        while line1:
            total_line += 1
            line1 = f1.readline()

            tmp_lig = line1[0:10]

            if tmp_lig in ligation_barcode_list:
                ligation_bc_match = ligation_barcode_list[tmp_lig]
                if ligation_bc_match in LIGwells_dict:
                    ligindx = LIGwells_dict[ligation_bc_match]
                    target_RT = line1[len(ligation_bc_match) + 14: len(ligation_bc_match) + 24]

                    if target_RT in RT_barcode_list:
                        barcode = RT_barcode_list[target_RT]
                        if barcode in RTsamplename_dict:
                            samplename = RTsamplename_dict[barcode]
                            RTwell = RTwells_dict[barcode]
                            filtered_line += 1
                            UMI = line1[len(ligation_bc_match) + 6: len(ligation_bc_match) + 14]
                            first_line = f'@LIG-{ligindx}_RT-{RTwell},PCR-{pcrwell},{UMI},{samplename}\n'
                            f3.write(first_line)

                            second_line = f2.readline()
                            f3.write(second_line)

                            third_line = f2.readline()
                            f3.write(third_line)

                            four_line = f2.readline()
                            f3.write(four_line)

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
            else:
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()

            line1 = f1.readline()
            line1 = f1.readline()
            line1 = f1.readline()

    print(f"PCR well: {pcrwell}, total line: {total_line}, filtered line: {filtered_line}, filter rate: {float(filtered_line) / float(total_line)}")


def attach_UMI_files(input_folder, pcrwell, output_folder, ligation_barcode_file, LIGlist, RT_barcode_file, RT_sample, core):
    init_message = f'''
    --------------------------start attaching UMI-----------------------------
    input folder: {input_folder}
    PCR well list: {pcrwell}
    output_folder: {output_folder}
    ligation barcode file: {ligation_barcode_file}
    ligation list: {LIGlist}
    RT barcode file: {RT_barcode_file}
    RT samplesheet: {RT_sample}
    ______________________________________________________________________________________________________________________________________________________
    '''
    print(init_message)

    print("Load ligation barcode dictionary...")

    with open(ligation_barcode_file, "rb") as barcodes:
        ligation_barcode_list = pickle.load(barcodes)

    print("Load RT barcode dictionary...")

    with open(RT_barcode_file, "rb") as barcodes:
        RT_barcode_list = pickle.load(barcodes)

    print("Loading sample names...")
    with open(RT_sample, "r") as RTsamplenames:
        RTreader = csv.DictReader(RTsamplenames)
        RTsamplename_dict = {rows["RTindex"]: rows["SampleName"] for rows in RTreader}
        samplelist = list(set(val for dic in RTsamplename_dict for val in RTsamplename_dict.values()))
        print(samplelist)

    with open(RT_sample, "r") as RTwells:
        RTreader3 = csv.DictReader(RTwells)
        RTwells_dict = {rows["RTindex"]: rows["RTwell"] for rows in RTreader3}
        RTwellslist = list(set(val for dic in RTwells_dict for val in RTwells_dict.values()))
        print(RTwells_dict)

    with open(LIGlist, "r") as LIGwells:
        RTreader4 = csv.DictReader(LIGwells)
        LIGwells_dict = {rows["LIGindex"]: rows["LIGwell"] for rows in RTreader4}

    pcrwell_file = open(pcrwell)
    pcrwell_list = []
    for line in pcrwell_file:
        pcrwell = line.strip()
        pcrwell_list.append(pcrwell)
    pcrwell_file.close()

    p = Pool(processes=int(core))
    func = partial(
        UMI_attach_read2_barcode_list, input_folder=input_folder, output_folder=output_folder, ligation_barcode_list=ligation_barcode_list,
        RT_barcode_list=RT_barcode_list, RTsamplename_dict=RTsamplename_dict, RTwells_dict=RTwells_dict, LIGwells_dict=LIGwells_dict, mismatch_rate=1
    )
    result = p.map(func, pcrwell_list)
    p.close()
    p.join()

    com_message = '''~~~~~~~~~~~~~~~UMI attachment done~~~~~~~~~~~~~~~~~~'''
    print(com_message)


if __name__ == "__main__":
    input_folder = sys.argv[1]
    pcrwell = sys.argv[2]
    output_folder = sys.argv[3]
    ligation_barcode_file = sys.argv[4]
    LIGlist = sys.argv[5]
    RT_barcode_file = sys.argv[6]
    RT_sample = sys.argv[7]
    core = sys.argv[8]
    attach_UMI_files(input_folder, pcrwell, output_folder, ligation_barcode_file, LIGlist, RT_barcode_file, RT_sample, core)
    
