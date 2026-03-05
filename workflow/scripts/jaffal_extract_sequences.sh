#!/usr/bin/env bash
# Extract fusion candidate sequences from raw FASTA based on gene counts
# Usage: Called by Snakemake with appropriate input/output parameters

set -euo pipefail

# Redirect all output to log file
exec > "${snakemake_log[0]}" 2>&1
echo "$(date)"

# Extract gene IDs from gene count file
cat "${snakemake_input[gene_count]}" | awk '{print $1}' > "${snakemake_output[fusion_temp]}"

# Extract sequences for candidate genes
reformat.sh "in=${snakemake_input[raw_fasta]}" out=stdout.fasta fastawrap=0 | \
    extract_seq_from_fasta "${snakemake_output[fusion_temp]}" > "${snakemake_output[fusion_fa]}"
