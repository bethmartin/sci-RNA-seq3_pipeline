#!/usr/bin/bash -l

#added another A so it wouldn't trim my synthetic A-tail 2022-03-07

input_folder=$1
sample=$2
output_folder=$3


echo Trimming sample: $sample
trim_galore $input_folder/$sample*.gz -a AAAAAAAAA --three_prime_clip_R1 1 -o $output_folder
echo Trimming $sample done.
