
input_folder=$1
sample=$2
output_folder=$3
barcode_file=$4
cutoff=$5

mismatch=1
python="/net/shendure/vol12/projects/sciRNAseq_script/anaconda2/bin/python2.7"
python_script="/net/shendure/vol1/home/martin91/scripts/sciRNAseq3/sam_split.py"

$python $python_script $input_folder/$sample.sam $barcode_file $output_folder $cutoff
echo splitting sample done: $sample
