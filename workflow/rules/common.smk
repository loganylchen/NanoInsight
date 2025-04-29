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
    if os.path.exists(f"data/{wildcards.sample}/fastq/pass.fq.gz"):
        return f"data/{wildcards.sample}/fastq/pass.fq.gz"
    else:
        return f"data/{wildcards.sample}/fastq/pass.fastq.gz"


def get_raw_blow5(wildcards):
    if os.path.exists(f"data/{wildcards.sample}/blow5/nanopore.blow5"):
        return f"data/{wildcards.sample}/blow5/nanopore.blow5"
    else:
        return f"data/{wildcards.sample}/blow5/nanopore.drs.blow5"





def get_final_output():
    # detect fusions
    final_output = []
    final_output += expand(
                '{project}/{sample}/alignment/{sample}_vg_genome_4aeron.gam', sample=list(samples.index),project=PROJECT)
    for tool in config['fusion']:
        if config['fusion'][tool]:
            final_output += expand(
                '{PROJECT}/{sample}/fusion/{tool}/{sample}_fusion_{tool}.tsv', sample=list(samples.index),PROJECT=PROJECT,tool=tool)
    return final_output