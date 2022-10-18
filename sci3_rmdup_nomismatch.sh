
input_folder=$1
sample=$2
output_folder=$3
mismatch=$4

python="/net/shendure/vol12/projects/sciRNAseq_script/anaconda2/bin/python2.7"
python_script="/net/shendure/vol1/home/martin91/scripts/sciRNAseq3/rm_dup_barcode_UMI_no_mismatch.py"

echo Filtering sample: $sample

$python $python_script $input_folder/$sample.sam $output_folder/$sample.sam $mismatch

echo Filtering $sample done.
