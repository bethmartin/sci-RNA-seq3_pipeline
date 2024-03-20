
#in this script, I will read into the sam file and the barcode file, and then count the reads number per barcode

#converted to python3 by ChatGPT :)

import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

def samfile_barcode_count(sam_file, barcode_file):
    # Generate the barcode list and barcode dictionary
    with open(barcode_file) as barcodes:
        barcode_ls = []
        barcode_dic = {}
        for line in barcodes:
            barcode = line.strip()
            barcode_ls.append(barcode)
            barcode_dic[barcode] = 0

    # Read the sam file, and count the number per barcode
    with open(sam_file) as sam:
        for line in sam:
            if line.startswith('@'):
                continue
            else:
                name = ((line.split('\t'))[0]).split(',')
                barcode = name[0]
                barcode_dic[barcode] += 1

    return barcode_dic

def split_samfile(sam_file, barcode_file, output_folder, cutoff):
    '''
    This script accepts a sam file, a barcode file, an output_file, a cutoff value,
    then it will call the samfile_barcode_count function and get the total read count per barcode,
    then it uses the cutoff value to filter the barcode,
    and generates the output samfile for single cells, generates the sample_ID.txt in the output folder,
    generates the reads distribution in the output folder/read_distribution_barcode;
    '''

    # Generate the count per barcode
    barcode_count = samfile_barcode_count(sam_file, barcode_file)

    # Plot the barcode reads distribution and save the result to the output folder
    plot_name = (sam_file.split('/')[-1]).split('.')[0]
    fig, ax = plt.subplots()
    ax.hist(barcode_count.values(), bins=100)
    ax.set_ylabel('frequency')
    ax.set_xlabel('Number of unique reads')
    fig_output = output_folder + '/' + plot_name + '.png'
    fig.savefig(fig_output)

    # Also output the barcode number and distribution to the output folder
    with open(output_folder + '/' + plot_name + '.txt', 'w') as read_dist:
        for barcode in barcode_count:
            line = barcode + ', %d\n' % (barcode_count[barcode])
            read_dist.write(line)

    # Filter the barcode based on the cutoff value
    barcode_filtered = [barcode for barcode in barcode_count if barcode_count[barcode] >= cutoff]

    # Generate the output sam file and sample_list file
    with open(output_folder + '/' + plot_name + '.' + 'sample_list.txt', 'w') as sample_list_file:
        output_files = {}
        for barcode in barcode_filtered:
            output_file = output_folder + '/' + plot_name + '.' + barcode + '.sam'
            output_files[barcode] = open(output_file, 'w')
            sample_list_file.write(plot_name + '.' + barcode + '\n')

        # Output each read to the output sam file
        with open(sam_file) as sam:
            for line in sam:
                if line.startswith('@'):
                    for barcode in barcode_filtered:
                        output_files[barcode].write(line)
                else:
                    barcode = ((line.split('\t'))[0]).split(',')[0]
                    if barcode in barcode_filtered:
                        output_files[barcode].write(line)

    # Close the files
    for barcode in barcode_filtered:
        output_files[barcode].close()

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("Usage: python script.py sam_file barcode_file output_folder cutoff")
        sys.exit(1)

    sam_file = sys.argv[1]
    barcode_file = sys.argv[2]
    output_folder = sys.argv[3]
    cutoff = int(sys.argv[4])
    split_samfile(sam_file, barcode_file, output_folder, cutoff)
    
