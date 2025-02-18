rule genion_fusion:
    input:
        fastq=get_raw_fastq,
        mapping_paf='{project}/{sample}/alignment/{sample}_minimap2_genome_4genion.paf',
        transcriptome_dup='{project}/resources/transcriptome_dup.tsv',
        dups='{project}/resources/genomicSuperDups.txt'
    output:
        fusions = '{project}/{sample}/fusion/genion/{sample}_fusion_genion.tsv'
    log:
        'logs/{project}/genion_fusion/{sample}.log'
    benchmark:
        'benchmarks/{project}/genion_fusion/{sample}.txt'
    params:
        gtf=config['reference']['annotation']
    resources:
        mem_mb = 1024 * 10
    threads: 1
    container:
        "docker://btrspg/genion:1.1.1"
    shell:
        "genion -i {input.fastq} "
        "-d {input.dups} "
        "--gtf {params.gtf} "
        "-g {input.mapping_paf} "
        "-s {input.transcriptome_dup} "
        "-o {output.fusions} 2>{log}"


