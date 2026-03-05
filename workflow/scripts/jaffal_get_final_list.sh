#!/usr/bin/env bash
# Generate final fusion candidate list by combining genome alignment and transcript data
# Usage: Called by Snakemake with appropriate input/output parameters

set -euo pipefail

# Redirect all output to log file
exec > "${snakemake_log[0]}" 2>&1
echo "$(date)"

# Check if genome PSL file has content
if [ ! -s "${snakemake_input[genome_psl]}" ]; then
    # Create empty output files if no genome alignment
    touch "${snakemake_output[final_list]}"
else
    # Process genome alignment with JAFFA R scripts
    R --vanilla --args \
        "${snakemake_input[genome_psl]}" \
        "${snakemake_input[fasta_reads_table]}" \
        "${snakemake_params[transtable]}" \
        "${snakemake_params[knowntable]}" \
        10000 \
        "NoSupport,PotentialReadThrough" \
        50 \
        "${snakemake_output[final_list]}" < /opt/softwares/JAFFA/make_final_table.R
fi

# Generate 3-gene fusion summary table
make_3_gene_fusion_table \
    "${snakemake_output[final_list]}" \
    "${snakemake_input[gene_count]}" \
    "${snakemake_output[three_gene_reads]}" > "${snakemake_output[three_gene_summary]}"
