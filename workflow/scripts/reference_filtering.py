#!/usr/bin/env python3
"""Reference file filtering utilities for NanoInsight workflow.

This module provides functions to filter FASTA and GTF reference files
to include only specified contigs/chromosomes.
"""

from typing import Set
from pathlib import Path
from Bio import SeqIO


def filter_fasta(
    input_fasta: str,
    output_fasta: str,
    contigs: Set[str]
) -> None:
    """Filter a FASTA file to keep only specified contigs.
    
    Args:
        input_fasta: Path to the input FASTA file.
        output_fasta: Path where the filtered FASTA will be written.
        contigs: Set of contig IDs to keep in the output.
        
    Raises:
        FileNotFoundError: If the input FASTA file doesn't exist.
        PermissionError: If unable to write to output path.
        
    Example:
        >>> filter_fasta(
        ...     "genome.fa",
        ...     "filtered_genome.fa",
        ...     {"1", "2", "X"}
        ... )
    """
    input_path = Path(input_fasta)
    if not input_path.exists():
        raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")
    
    with open(input_fasta) as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id in contigs:
                SeqIO.write(record, outfile, "fasta")


def filter_gtf(
    input_gtf: str,
    output_gtf: str,
    contigs: Set[str]
) -> None:
    """Filter a GTF file to keep only entries from specified contigs.
    
    Args:
        input_gtf: Path to the input GTF file.
        output_gtf: Path where the filtered GTF will be written.
        contigs: Set of contig IDs to keep in the output.
        
    Raises:
        FileNotFoundError: If the input GTF file doesn't exist.
        PermissionError: If unable to write to output path.
        
    Note:
        Comment lines (starting with '#') are always preserved.
    """
    input_path = Path(input_gtf)
    if not input_path.exists():
        raise FileNotFoundError(f"Input GTF file not found: {input_gtf}")
    
    with open(input_gtf) as infile, open(output_gtf, "w") as outfile:
        for line in infile:
            if line.startswith("#") or line.split("\t")[0] in contigs:
                outfile.write(line)


# Snakemake integration - these variables are injected by Snakemake
if __name__ == "__main__":
    import snakemake
    
    input_fasta = snakemake.input.fasta
    input_gtf = snakemake.input.gtf
    output_fasta = snakemake.output.fasta
    output_gtf = snakemake.output.gtf
    contigs = set(snakemake.params.contigs)
    
    filter_fasta(input_fasta, output_fasta, contigs)
    filter_gtf(input_gtf, output_gtf, contigs)
