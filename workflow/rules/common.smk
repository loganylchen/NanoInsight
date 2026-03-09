import glob
import pandas as pd
import sys
import os
from snakemake.utils import validate
from snakemake.logging import logger

PROJECT = config["project"]

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"SampleName": str}, comment="#")
    .set_index("SampleName", drop=False)
    .sort_index()
)



def get_raw_fastq(wildcards):
    """Get the raw FASTQ file for a sample.
    
    Checks for both .fq.gz and .fastq.gz extensions.
    Raises FileNotFoundError if neither exists.
    """
    fq_gz = f"data/{wildcards.sample}/fastq/pass.fq.gz"
    fastq_gz = f"data/{wildcards.sample}/fastq/pass.fastq.gz"
    
    if os.path.exists(fq_gz):
        return fq_gz
    elif os.path.exists(fastq_gz):
        return fastq_gz
    else:
        raise FileNotFoundError(
            f"No FASTQ file found for sample '{wildcards.sample}'. "
            f"Expected one of: {fq_gz} or {fastq_gz}"
        )




def get_raw_blow5(wildcards):
    """Get the raw BLOW5 file for a sample.
    
    Checks for both nanopore.blow5 and nanopore.drs.blow5 naming conventions.
    Raises FileNotFoundError if neither exists.
    """
    blow5 = f"data/{wildcards.sample}/blow5/nanopore.blow5"
    drs_blow5 = f"data/{wildcards.sample}/blow5/nanopore.drs.blow5"
    
    if os.path.exists(blow5):
        return blow5
    elif os.path.exists(drs_blow5):
        return drs_blow5
    else:
        raise FileNotFoundError(
            f"No BLOW5 file found for sample '{wildcards.sample}'. "
            f"Expected one of: {blow5} or {drs_blow5}"
        )





def get_final_output():
    # detect fusions
    final_output = []
    final_output += expand(
                "{project}/{sample}/fusion/aeron/{sample}_tofusion.sam", sample=list(samples.index),project=PROJECT)
    
    for tool in config['fusion']:
        if config['fusion'][tool]:
            final_output += expand(
                '{PROJECT}/{sample}/fusion/{tool}/{sample}_fusion_{tool}.tsv', sample=list(samples.index),PROJECT=PROJECT,tool=tool)
    return final_output