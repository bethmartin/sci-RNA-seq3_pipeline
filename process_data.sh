# Updated on Apr-6-2024, by CX Qiu

#!/bin/bash
# this scRNA-seq pipeline accept a input folder, and then use the default parameter for the data processing and analysis

run_ID="<put your run ID here>"

my_folder="<put your working folder here>"

sci_experiment_type="MEGA"

num=${SGE_TASK_ID}

# define the fastq folder including all fastq files
fastq_folder="${my_folder}/${run_ID}/nobackup/fastq"

#make sure not to put a slash on the end of this path
work_folder="${my_folder}/${run_ID}/nobackup/output_${num}"

# define the PCR group id after demultiplexing
sample_ID="${my_folder}/${run_ID}/sample_ID_${num}.txt"

RT_sample="${my_folder}/${run_ID}/RTsamplesheet.csv"

# define the core number for parallel processing
core=12
core_sam=3

# define the number of unique reads cutoff for splitting single cell
cutoff=100

# mouse genome
index="/net/shendure/vol10/nobackup/genome/STAR/GRCm39-p6-all-PEmax"
gtf_file="/net/shendure/vol10/nobackup/genome/GTF/gencode.vM37.chr_patch_hapl_scaff.annotation_add_PEmax.gtf.gz"

#define the mismatch rate for removing duplicates:
mismatch=1

#define the bin of python (python V2.7)
python_path="/net/shendure/vol12/projects/sciRNAseq_script/anaconda2/bin/"

#define the bin of R (R V3.6.3)
R_path="/net/gs/vol3/software/modules-sw/R/4.3.2/Linux/Ubuntu22.04/x86_64/bin"

#define the location of script:
script_path="/net/gs/vol1/home/cxqiu/work/scripts/JAX_rna/Jun_pipeline"
# define the location of the R script for multi-core processing
Jun_R_script=$script_path/sci3_bash_input_ID_output_core.R

if [ "$sci_experiment_type" = "MEGA" ]; then
    # define the location of the ligation barcodes
    ligation_barcode=$script_path/lig_MEGA_bc.pickle2
    # define the location of the RT barcodes
    RT_barcode=$script_path/RT_MEGA_bc.pickle2
    # define the location of the combined RT and ligation barcodes
    barcodes=$script_path/combined_768_wells.txt
    ligwells=$script_path/ligationwellbarcodeMEGA.csv
else
    # define the location of the ligation barcodes
    ligation_barcode=$script_path/lig_384_bc.pickle2
    # define the location of the RT barcodes
    RT_barcode=$script_path/RT_384_bc.pickle2
    # define the location of the combined RT and ligation barcodes
    barcodes=$script_path/combined_384_bc.txt
fi

now=$(date)
echo "Current time : $now"

module load modules modules-init modules-gs
module load samtools/1.19
module load bedtools/2.31.1

UMI_attach_folder=$work_folder/UMI_attach
trimmed_fastq_folder=$work_folder/trimmed_fastq
STAR_alignment_folder=$work_folder/STAR_alignment
filtered_sam_folder=$work_folder/filtered_sam
rmdup_sam_folder=$work_folder/rmdup_sam
sam_splitted_folder=$work_folder/sam_splitted
report_folder=$work_folder/report

##########################
### Step 1: UMI attach ###
##########################
### ignoring I1 and I2 (PCR well information)
### adding R1 (barcode information) to R2 (mRNA)

echo "changing the name of the fastq files..."
for sample in $(cat $sample_ID); do
    echo changing name $sample
    mv $fastq_folder/*$sample*R1*.gz $fastq_folder/$sample.R1.fastq.gz
    mv $fastq_folder/*$sample*R2*.gz $fastq_folder/$sample.R2.fastq.gz
done

echo "Attaching barcode and UMI...."
mkdir -p $UMI_attach_folder

if [ "$sci_experiment_type" = "MEGA" ]; then
    $python_path/python \
    $script_path/UMI_barcode_attach_gzipped_with_dic_MEGA.py \
    $fastq_folder \
    $sample_ID \
    $UMI_attach_folder \
    $ligation_barcode \
    $ligwells \
    $RT_barcode \
    $RT_sample \
    $core
else
    $python_path/python \
    $script_path/UMI_barcode_attach_gzipped_with_dic.py \
    $fastq_folder \
    $sample_ID \
    $UMI_attach_folder \
    $ligation_barcode \
    $RT_barcode \
    $core
fi

### calculate read number in the UMI_attach folder
for i in `ls $UMI_attach_folder/*.R2.fastq.gz`; do
    echo $(( $(zcat "$i" | wc -l) / 4 )) >> $work_folder/read_num_UMI_attach.txt
done

echo ">>>Step 1 (UMI attach) has been done."


##################################
### Step 2: Trimming the read2 ###
##################################

echo
echo "Start trimming the read2 file..."

module load python/3.12.1
module load cutadapt/4.6
module load fastqc/0.12.1
module load trim_galore/0.6.10

$R_path/Rscript \
    $Jun_R_script \
    $script_path/sci3_trim_beth.sh \
    $UMI_attach_folder \
    $sample_ID \
    $trimmed_fastq_folder \
    $core

module unload python/3.12.1

echo ">>>Step 2 (trmming read2 file) has been done."


##############################
### Step 3: Aligning reads ###
##############################

echo
echo "Start aligning reads..."

module load STAR/2.6.1d

mkdir -p $STAR_alignment_folder
STAR --genomeDir $index --genomeLoad Remove
for sample in $(cat $sample_ID); do 
    echo Aligning $sample
    STAR \
    --runThreadN $core \
    --outSAMstrandField intronMotif \
    --genomeDir $index \
    --readFilesCommand zcat \
    --readFilesIn $trimmed_fastq_folder/$sample*gz \
    --outFileNamePrefix $STAR_alignment_folder/$sample \
    --genomeLoad LoadAndKeep \
    --outReadsUnmapped Fastx 
done
STAR --genomeDir $index --genomeLoad Remove

echo ">>>Step 3 (aligning reads) has been done."


#############################
### Step 4: Filtering SAM ###
#############################

echo
echo "Start filter and sort the sam files..."

$R_path/Rscript \
    $Jun_R_script \
    $script_path/sci3_filter.sh \
    $STAR_alignment_folder \
    $sample_ID \
    $filtered_sam_folder \
    $core_sam

echo ">>>Step 4 (filtering sam) has been done."

### deleted the UMI_attach, trimmed reads, and STAR output
### rm -rf $UMI_attach_folder/
rm -rf $trimmed_fastq_folder/
rm -rf $STAR_alignment_folder/


###################################
### Step 5: removing duplicates ###
###################################

echo
echo "Start removing duplicates..."
mkdir -p $rmdup_sam_folder
module unload python

$R_path/Rscript \
    $Jun_R_script \
    $script_path/sci3_rmdup_nomismatch_beth.sh \
    $filtered_sam_folder \
    $sample_ID \
    $rmdup_sam_folder \
    $core \
    $mismatch

mkdir -p $report_folder/duplicate_read
mv $rmdup_sam_folder/*.csv $report_folder/duplicate_read/

### calculate read number and estimate duplicate rate
for i in `ls $filtered_sam_folder/*.sam`; do 
    samtools view -c $i >> $work_folder/tmp1
done
for i in `ls $rmdup_sam_folder/*.sam`; do 
    samtools view -c $i >> $work_folder/tmp2
done
paste $work_folder/tmp1 $work_folder/tmp2 > $work_folder/read_num.txt
rm $work_folder/tmp1 $work_folder/tmp2

echo ">>>Step 5 (removing duplicates) has been done."


##################################
### Step 6: split the sam file ###
##################################

echo
echo "Start splitting the sam file..."
module unload python

$R_path/Rscript \
    $Jun_R_script \
    $script_path/sci3_split_beth.sh \
    $rmdup_sam_folder \
    $sample_ID \
    $sam_splitted_folder \
    $core \
    $barcodes \
    $cutoff

cat $sam_splitted_folder/*sample_list.txt > $sam_splitted_folder/All_samples.txt
cp $sam_splitted_folder/All_samples.txt $work_folder/barcode_samples.txt

mkdir -p $report_folder/barcode_read_distribution
mv $sam_splitted_folder/*.txt $report_folder/barcode_read_distribution/
mv $sam_splitted_folder/*.png $report_folder/barcode_read_distribution/

echo ">>>Step 6 (splitting the sam file) has been done."


#################################
### Step 7: create gene count ###
#################################

gene_count_folder=$report_folder/human_mouse_gene_count
sample_ID=$work_folder/barcode_samples.txt

echo "Start the gene count...."
$python_path/python \
    $script_path/sciRNAseq_count.py \
    $gtf_file \
    $sam_splitted_folder \
    $sample_ID \
    $core

echo "Make the output folder and transfer the files..."
mkdir -p $gene_count_folder
find $sam_splitted_folder -name *.count -exec cat {} + > $gene_count_folder/count.MM
find $sam_splitted_folder -name '*.count' | xargs rm -f
find $sam_splitted_folder -name *.report -exec cat {} + > $gene_count_folder/report.MM
find $sam_splitted_folder -name '*.report' | xargs rm -f
mv $sam_splitted_folder/*_annotate.txt $gene_count_folder/

echo ">>>Step 7 (creating gene counts) has been done."


##################################
### Step 8: create data matrix ###
##################################

echo "Start creating data matrix...."
$R_path/Rscript \
    $script_path/gene_count_processing_sciRNAseq.R \
    $report_folder

echo ">>>Step 8 (creating data matrix) has been done."

### rm -rf $filtered_sam_folder/
rm -rf $rmdup_sam_folder/
rm -rf $sam_splitted_folder/

now=$(date)
echo "Current time : $now"
echo ">>> Hello, world!"

