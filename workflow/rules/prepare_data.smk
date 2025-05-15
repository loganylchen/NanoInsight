rule prep_raw_fastq:
    input:
        fastq=get_raw_fastq
    output:
        fastq=temp("{project}/{sample}/raw_fastq/{sample}_raw.fastq.gz")
    log:
        "logs/{project}/prep_raw_fastq/{sample}.log"
    threads: 1
    resources:
        mem_mb = 1024 * 5
    container:
        config['containers']['aeron']
    benchmark:
        "benchmarks/{project}/prep_raw_fastq/{sample}.txt"
    script:
        '../scripts/prep_raw_fastq.sh'


rule prep_unzip_fastq:
    input:
        fastq="{project}/{sample}/raw_fastq/{sample}_raw.fastq.gz"
    output:
        fastq="{project}/{sample}/raw_fastq/{sample}_raw.fastq"
    log:
        "logs/{project}/prep_unzip_fastq/{sample}.log"
    threads: 1
    resources:
        mem_mb = 1024 * 5
    container:
        config['containers']['aeron']
    benchmark:
        "benchmarks/{project}/prep_unzip_fastq/{sample}.txt"
    shell:
        'gzip -dc {input.fastq} > {output.fastq} && echo `date` > {log} '


