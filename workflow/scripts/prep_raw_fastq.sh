#!/usr/bin/env bash
# Preprocess raw FASTQ files for the NanoInsight pipeline.
# Removes extra fields from read headers (keeps only the read ID).
#
# Input: Raw FASTQ file (gzipped)
# Output: Cleaned FASTQ file with simplified headers (gzipped)

set -euo pipefail

exec > "${snakemake_log[0]}" 2>&1
echo "Starting prep_raw_fastq at $(date)"

input_fastq="${snakemake_input[fastq]}"
output_fastq="${snakemake_output[fastq]}"

echo "Input: ${input_fastq}"
echo "Output: ${output_fastq}"

# Process FASTQ: keep only read ID (first field) from headers
zcat "${input_fastq}" | awk -F ' ' '{print $1}' | gzip -c > "${output_fastq}"

echo "Input line count:"
zcat "${input_fastq}" | wc -l

echo "Output line count:"
zcat "${output_fastq}" | wc -l

echo "Completed at $(date)"
