#!/usr/bin/env bash

set -x
set -e


exec > "${snakemake_log[0]}" 2>&1
echo `date` 

if [ ! -s ${snakemake_input[genome_psl]} ]; then
    touch ${snakemake_output[final_list]};
else
    R --vanilla --args ${snakemake_input[genome_psl]} \
        ${snakemake_input[fasta_reads_table]} ${snakemake_params[transtable]} \
        ${snakemake_params[knowntable]} 10000 \
        "NoSupport,PotentialReadThrough" 50 \
        ${snakemake_output[final_list]} < /opt/softwares/JAFFA/make_final_table.R
fi

make_3_gene_fusion_table ${snakemake_output[final_list]} ${snakemake_input[gene_count]} ${snakemake_output[three_gene_reads]}  > ${snakemake_output[three_gene_summary]} 



