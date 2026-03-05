#!/usr/bin/env bash
# Compile final JAFFAL fusion results and generate output FASTA
# Usage: Called by Snakemake with appropriate input/output parameters

set -euo pipefail

# Redirect all output to log file
exec > "${snakemake_log[0]}" 2>&1
echo "$(date)"

# Run JAFFA compile results R script
R --vanilla --args \
    "${snakemake_params[final_out_prefix]}" \
    "${snakemake_input[three_gene_summary]}" \
    "${snakemake_input[final_list]}" < /opt/softwares/JAFFA/compile_results.R

# Extract fusion sequences for each candidate
while read -r line; do
    /opt/softwares/JAFFA/scripts/get_fusion_seqs.bash "$line" "${snakemake_output[output_fasta]}"
done < "${snakemake_params[final_out_prefix]}.csv"

echo "Done writing fusion sequences"
echo "All Done."
echo "*************************************************************************"
echo " Citation for JAFFA_direct, JAFFA_assembly and JAFFA_hybrid: "
echo "   Davidson, N.M., Majewski, I.J. & Oshlack, A. "
echo "   JAFFA: High sensitivity transcriptome-focused fusion gene detection."
echo "   Genome Med 7, 43 (2015)"
echo "*************************************************************************"
echo " Citation for JAFFAL: "
echo "   Davidson, N.M. et al. "
echo "   JAFFAL: detecting fusion genes with long-read transcriptome sequencing"
echo "   Genome Biol. 23, 10 (2022)"
echo "*************************************************************************"

# Create final output marker
touch "${snakemake_output[fusion]}"
