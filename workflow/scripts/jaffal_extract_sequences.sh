#!/usr/bin/env bash


set -x
set -e


exec > "${snakemake_log[0]}" 2>&1
echo `date` 


cat ${snakemake_input[gene_count]} | awk '{print $1}' > ${snakemake_output[fusion_temp]} ;
reformat.sh in=${snakemake_input[raw_fasta]} out=stdout.fasta fastawrap=0 | extract_seq_from_fasta ${snakemake_output[fusion_temp]} > ${snakemake_output[fusion_fa]}





