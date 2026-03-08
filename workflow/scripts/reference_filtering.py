#!/usr/bin/env python3
"""Filter genome FASTA and GTF files to include only specified contigs.

This script filters both a reference genome FASTA file and its corresponding
GTF annotation file to include only the chromosomes/contigs specified in the
configuration. This is useful for reducing memory usage and processing time
when only certain chromosomes are needed for analysis.

Usage:
    This script is designed to be called from Snakemake with the following
    inputs/outputs/params defined in the rule.
"""

import sys
from typing import List
from Bio import SeqIO


def filter_fasta(input_fasta: str, output_fasta: str, contigs: List[str]) -> None:
    """Filter FASTA file to include only specified contigs.
    
    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output filtered FASTA file
        contigs: List of contig/chromosome names to keep
    """
    with open(input_fasta) as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id in contigs:
                SeqIO.write(record, outfile, "fasta")


def filter_gtf(input_gtf: str, output_gtf: str, contigs: List[str]) -> None:
    """Filter GTF file to include only specified contigs.
    
    Args:
        input_gtf: Path to input GTF file
        output_gtf: Path to output filtered GTF file  
        contigs: List of contig/chromosome names to keep
    """
    with open(input_gtf) as infile, open(output_gtf, "w") as outfile:
        for line in infile:
            if line.startswith("#") or line.split("\t")[0] in contigs:
                outfile.write(line)


# Main execution when called from Snakemake
input_fasta = snakemake.input.fasta
input_gtf = snakemake.input.gtf
output_fasta = snakemake.output.fasta
output_gtf = snakemake.output.gtf
contigs = snakemake.params.contigs

filter_fasta(input_fasta, output_fasta, contigs)
filter_gtf(input_gtf, output_gtf, contigs)
