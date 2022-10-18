# SciRNAseq3 pipeline

originally written by Junyue Cao, tweaked by Beth

Modules used throughout pipeline:\
module load modules modules-init modules-gs\
module load gmp/6.1.2\
module load mpfr/4.0.1\
module load mpc/1.1.0\
module load gcc/8.2.0\
module load bcl2fastq/2.20\
module load STAR/2.6.1d\
module load samtools/1.9\
module load bedtools/2.27.1\
module load R/4.0.0\
python 2.7 installed somewhere

Several files will need editing to find your folders where your data and scripts are kept, when you are downloading these scripts for the first time.

This pipeline is typically run with 10 cores.

Before beginning the pipeline you will need to create a STAR reference to align the reads to, if you don't have one already.
The typical command for making this is like:

`STAR --runMode genomeGenerate --genomeDir <where you want your genome, must exist> --genomeFastaFiles <fasta1.fa fasta2.fa etc.> --sjdbGTFfile <gtf file that goes with the sequences, must be unzipped> --runThreadN 10`

Make sure you have requested enough memory/cores to hold the genome. 

### Demultiplexing the sequencing run
script is ***demultiplex.sh***\
In this script you will need to change the path to the run folder, samplesheet and the output folder. \
You will also change the --use-bases-mask argument depending on whether its a single- or double-indexed run\
A couple of samplesheet templates that include the 4 plates of PCR P7 indexes. (add your p5 index to the dual index version).\
***demux-samplesheet-singleindex.csv***\
***demux-samplesheet-dualindex.csv***

### Next is running the main sci-RNA-seq3 pipeline. Main script is scRNA_seq_pipeline.sh
If you've just downloaded the files, here are the changes you need to make:

***scRNA_seq_pipeline.sh*** \
    Change the "script_path" at the top to the folder where you've put all the scripts\
    Change the "python_path" in the common settings to where you have python v2.7\
    For *each* experiment, you will have to change the experiment-specific settings\
        fastq_folder = where you put the fastq files after demuxing\
        all_output_folder = where all the output is going\
        pcrwell = this is the same as the first field of your samplesheet (sample_ID), just listed one per line in a text file\
        RT_sample = a samplesheet with the headers: RTwell, RTindex, SampleName. Defines which samples went into each well of the RT plate(s).\
        index = the STAR reference folder\
        gtf_file = the .gtf file that goes with that reference (is gzipped)\
        countscript = sciRNAseq_count.py should work for most references
        

***sci3_rmdup.sh***\
    Change the folders for "python" and "python_script"

***sci3_rmdup_nomismatch.sh***\
    Change the folders for "python" and "python_script"
    
***sci3_split.sh***\
    Change the folders for "python" and "python_script"


### What's happening in the script

The very first thing it does is simplify the names of the read files. 

The next step is to take those read files and pull out the Ligation index, RT index, and UMI from the read1, and put this information into the @ line.\
This is tricky because the ligation index can either be 9 or 10bp long. (This is done to shift half the reads and add base complexity for the constant sequence downstream), so the script figures out which ligation index it's supposed to be first, so once it knows the length it will know the register to pull out the indexes and UMI.\
The script originally put the index sequences themselves into the @ line, but I've tweaked it to put the well addresses that the indexes come from, 
so that it is more human-readable. I've also added the sample name to the @ line as well, taking that info from the ***RTsamplesheet.csv***.

Then the reads are trimmed of their poly A tails. It looks for runs of A longer than 9bp on the 3' end.

The reads are then aligned with STAR

Reads are filtered for quality and sorted

Duplicates are removed based on exact UMI sequences, then again for 1 mismatch (this second filtering can probably be removed)

The sam files are split into files for each cell. This is when things get a little crazy with the file numbers, especially with large runs. 

Then the genes are counted. 

That's the end of the main script and then the final step is to do the gene count processing.\
    ***gene_count_processing_sciRNAseq_CX.R***\
        In this script, for each experiment, change the report folder (should be <your experiment folder>/nobackup/output/report) and output folder (this can be the same folder)\
        load R/4.0.0 and then run the genecount processing script with: `Rscript gene_count_processing_sciRNAseq_CX.R` 


### Extras

***estimate_dup_rate.sh***\
    This will compare the sizes of the sam files before and after filtering to help estimate your duplication rate. Run the script from your experiment 
    directory, then open up the resulting file in R with these lines to get your dup rate:\
    dat = read.table("read_num.txt", as.is=T)\
    dup = 1 - (dat$V3/dat$V1)\
    summary(dup)

***generate_RT_lig_pickle.py***
    you will only have to generate new pickle files if you add new custom indexes to the ones already used. \
    MAKING PICKLE DICTIONARY\
    script:\
    ***generate_RT_lig_pickle.py***\
    input files:\
    ***RT_gestalt_bc.txt***\
    ***ligation_gestalt_bc.txt***\
    output # change these paths in the script as needed\
    path_to_your_script_directory/lig_gestalt_bc.pickle\
    path_to_your_script_directory/lig_gestalt_bc.pickle2\
    path_to_your_script_directory/RT_gestalt_bc.pickle\
    path_to_your_script_directory/RT_gestalt_bc.pickle2\
    to run:\
    <path to python 2.7>/anaconda2/bin/python <path to this script>/generate_RT_lig_pickle.py    
