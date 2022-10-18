#added another A so it wouldn't trim my synthetic A-tail 2022-03-07

input_folder=$1
sample=$2
output_folder=$3

#had to tweak some modules for centos7, not sure if these will work, but trying out
module load python/3.7.7 
module load cutadapt/1.18  
module load trim_galore/0.6.5 
echo Trimming sample: $sample
trim_galore $input_folder/$sample*.gz -a AAAAAAAAA --three_prime_clip_R1 1 -o $output_folder
module unload python/3.7.7
echo Trimming $sample done.