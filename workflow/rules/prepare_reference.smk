rule extract_dup:
    output:
        '{project}/resources/genomicSuperDups.txt'
    params:
        dup_file_link=config['misc']['genion_dup_file_link']
    threads: 1
    resources:
        mem_mb = 1024 * 5
    priority: 1
    log:
        'logs/{project}/prepre_dups.log'
    shell:
        "wget -q -O - {params.dup_file_link} | zcat | tee {output} 1>{log}"



rule extract_transcripts:
    input:
        fasta="{project}/resources/genome.fasta",
        gtf="{project}/resources/genome.gtf"
    output:
        transcripts="{project}/resources/transcriptome.fasta"
    log:
        "logs/{project}/ref/get_transcripts.log"
    container:
        "docker://btrspg/gffread:0.12.7"
    benchmark:
        "benchmarks/{project}/extract_transcripts.benchmark.txt"
    shell:
        "gffread -w {output.transcripts} -g {input.fasta} {input.gtf} 2> {log}"



rule filtering_genome_and_annotation:
    input:
        fasta=config['reference']['genome'],
        gtf=config['reference']['annotation'],
    output:
        fasta="{project}/resources/genome.fasta",
        gtf="{project}/resources/genome.gtf",
    log:
        "logs/{project}/ref/filtering_references.log",
    params:
        contigs=config["reference"]["contigs"],
    container:
        "docker://btrspg/biopython:1.85"
    benchmark:
        "benchmarks/{project}/filtering_references.benchmark.txt"
    script:
        "../scripts/reference_filtering.py"



rule vg_genome_index_aeron:
    input:
        reference = "{project}/resources/genome.fasta",
        annotation= "{project}/resources/genome.gtf",
    output:
        '{project}/resources/genome.gfa',
    log:
        'logs/{project}/vg_genome_index_aeron.log'
    container:
        config['containers']['aeron']
    threads: 1
    benchmark:
        'benchmarks/{project}/vg_genome_index_aeron.txt'
    resources:
        mem_mb = 1024 * 300,
    shell:
        "GraphBuilder.py -e {input.reference} -g {input.annotation} -o {output} 2>{log} 1>&2  "