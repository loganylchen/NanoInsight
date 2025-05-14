#!/usr/bin/env bash




exec > "${snakemake_log[0]}" 2>&1
echo `date` 

set -e
input_fastq="${snakemake_input[fastq]}"
output_fastq="${snakemake_output[fastq]}"



zcat ${input_fastq} | awk -F ' ' '{print $1}' | gzip -c > ${output_fastq} 
zcat ${input_fastq} | wc -l 
zcat ${input_fastq} | head -n 4
echo `date`
zcat ${output_fastq}  | wc -l 
zcat ${output_fastq}  | head -n 4



