#!/usr/bin/env bash

set -x 
set -e

exec > "${snakemake_log[0]}" 2>&1
date

grep $'\t+\t' "${snakemake_input[0]}" | awk -F'\t' -v OFS="\t" '{ print $4-$3,0,0,0,0,0,0,0,$5,$1,$2,$3,$4,$6,$7,$8,$9,2, 100","$4-$3-100",",$3","$3+100",",  $8","$9-$4+$3+100"," }' > "${snakemake_output[0]}" ;
grep $'\t-\t' "${snakemake_input[0]}" | awk -F'\t' -v OFS="\t" '{ print $4-$3,0,0,0,0,0,0,0,$5,$1,$2,$3,$4,$6,$7,$8,$9,2, 100","$4-$3-100",", $2-$4","$2-$4+100",", $8","$9-$4+$3+100"," }' >> "${snakemake_output[0]}" ;