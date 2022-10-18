#from Ronnie
#!/bin/bash

#change these 3 things, make sure to put your full path, also change --use-bases-mask in the command line as needed for single vs dual indexed runs
run_folder="/net/shendure/vol9/seq/NEXTSEQ/??????"
sample_sheet="<your experiment folder>/demux-samplesheet-singleindex.csv"
output_folder="<your experiment folder>/nobackup/fastq"

module load modules modules-init modules-gs
module load gmp/6.1.2
module load mpfr/4.0.1
module load mpc/1.1.0
module load gcc/8.2.0
module load bcl2fastq/2.20

echo "---------------start demultiplex-----------------------"
echo $(date)
echo "run folder is $run_folder"
echo "sample sheet is $sample_sheet"
echo "output_folder is $output_folder"

#create the output folder:
echocr	
echo start making the output_folder
mkdir -p $output_folder/report


# do the demultiplex
#change bases mask as needed ie --use-bases-mask Y*,I*,Y* for singleindex, --use-bases-mask Y*,I*,I*,Y* for dual index
bcl2fastq --runfolder-dir $run_folder -o $output_folder --sample-sheet $sample_sheet --reports-dir $output_folder/report --barcode-mismatches 1 --create-fastq-for-index-reads --no-lane-splitting --use-bases-mask Y*,I*,Y* --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0

echo "------------------demultiplex done -----------------"