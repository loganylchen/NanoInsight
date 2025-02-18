rule longgf_fusion:
    input:
        mapping_bam='{project}/{sample}/alignment/{sample}_minimap2_genome_4longgf.bam',
    output:
        fusions = '{project}/{sample}/fusion/longgf/{sample}_fusion_longgf.tsv'
    params:
        gtf=config['reference']['annotation']
    log:
        'logs/{project}/longgf_fusion/{sample}.log'
    benchmark:
        'benchmarks/{project}/longgf_fusion/{sample}.txt'
    resources:
        mem_mb = 1024 * 10
    threads: 1
    container:
        "docker://btrspg/longgf:0.1.2"
    shell:
        'LongGF {input.mapping_bam}  {params.gtf} 100 30 100 '
        '| grep "SumGF" > {output.fusions}'


