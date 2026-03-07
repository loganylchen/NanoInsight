"""
Reference filtering script for NanoInsight workflow.

Filters FASTA and GTF files to keep only specified contigs.
Used to reduce reference size and focus analysis on canonical chromosomes.
"""

import sys
from typing import List
from Bio import SeqIO


def filter_fasta(input_fasta: str, output_fasta: str, contigs: List[str]) -> None:
    """
    Filter a FASTA file to include only specified contigs.

    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output filtered FASTA file
        contigs: List of contig names to keep
    """
    with open(input_fasta) as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id in contigs:
                SeqIO.write(record, outfile, "fasta")


def filter_gtf(input_gtf: str, output_gtf: str, contigs: List[str]) -> None:
    """
    Filter a GTF file to include only specified contigs.

    Args:
        input_gtf: Path to input GTF file
        output_gtf: Path to output filtered GTF file
        contigs: List of contig names to keep
    """
    with open(input_gtf) as infile, open(output_gtf, "w") as outfile:
        for line in infile:
            if line.startswith("#") or line.split("\t")[0] in contigs:
                outfile.write(line)


# Main execution when run as snakemake script
input_fasta: str = snakemake.input.fasta
input_gtf: str = snakemake.input.gtf
output_fasta: str = snakemake.output.fasta
output_gtf: str = snakemake.output.gtf
contigs: List[str] = snakemake.params.contigs

filter_fasta(input_fasta, output_fasta, contigs)
filter_gtf(input_gtf, output_gtf, contigs)
