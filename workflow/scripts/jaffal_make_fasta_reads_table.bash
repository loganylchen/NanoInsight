#!/usr/bin/env bash


set -x
set -e


exec > "${snakemake_log[0]}" 2>&1
echo `date` 


echo  -e "transcript\tbreak_min\tbreak_max\tfusion_genes\tspanning_pairs\tspanning_reads" > ${snakemake_output[fasta_reads_table]}
awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"0"\t"1}' ${snakemake_input[gene_count]} | sort -u  >> ${snakemake_output[fasta_reads_table]}


