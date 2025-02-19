rule minimap2_genome_fusion_genion:
    input:
        fastq=get_raw_fastq
    output:
        mapping_paf='{project}/{sample}/alignment/{sample}_minimap2_genome_4genion.paf'
    params:
        reference=config['reference']['genome'],
        ext_opt=config['params']['minimap2_genion'],
    threads: config['threads']['minimap2']
    resources:
        mem_mb = 1024 * 30
    log:
        'logs/{project}/minimap2_genome_fusion_genion/{sample}.log'
    benchmark:
        'benchmarks/{project}/minimap2_genome_fusion_genion/{sample}.txt'
    conda:
        '../envs/minimap2.yaml'
    shell:
        "minimap2 {params.ext_opt} -t {threads}  {params.reference}  {input.fastq} -o {output.mapping_paf} 2>{log}"


rule minimap2_genome_fusion_longgf:
    input:
        fastq=get_raw_fastq
    output:
        mapping_bam='{project}/{sample}/alignment/{sample}_minimap2_genome_4longgf.bam'
    params:
        reference=config['reference']['genome'],
        ext_opt=config['params']['minimap2_longgf'],
    threads: config['threads']['minimap2']
    resources:
        mem_mb = 1024 * 30
    log:
        'logs/{project}/minimap2_genome_fusion_longgf/{sample}.log'
    benchmark:
        'benchmarks/{project}/minimap2_genome_fusion_longgf/{sample}.txt'
    conda:
        '../envs/minimap2.yaml'
    shell:
        "minimap2 {params.ext_opt} -t {threads}  {params.reference}  {input.fastq}  2>{log} "
        "|samtools view -bSh "
        "| samtools sort -n - -o {output.mapping_bam}"


rule minimap2_transcriptome_fusion_jaffal:
    input:
        fastq=get_raw_fastq
    output:
        mapping_paf='{project}/{sample}/alignment/{sample}_minimap2_transcriptome_4jaffal.paf'
    params:
        reference=config['reference']['transcriptome'],
    threads: config['threads']['minimap2']
    resources:
        mem_mb = 1024 * 30
    log:
        'logs/{project}/minimap2_transcriptome_fusion_jaffal/{sample}.log'
    benchmark:
        'benchmarks/{project}/minimap2_transcriptome_fusion_jaffal/{sample}.txt'
    conda:
        '../envs/minimap2.yaml'
    shell:
        "minimap2 -x map-ont -t {threads} -c  {params.reference}  {input.fastq} -o {output.mapping_paf} 2>{log}"


rule minimap2_genome_fusion_jaffal:
    input:
        fastq=get_raw_fastq
    output:
        mapping_paf='{project}/{sample}/alignment/{sample}_minimap2_genome_4jaffal.paf'
    params:
        reference=config['reference']['genome'],
    threads: config['threads']['minimap2']
    resources:
        mem_mb = 1024 * 30
    log:
        'logs/{project}/minimap2_genome_fusion_jaffal/{sample}.log'
    benchmark:
        'benchmarks/{project}/minimap2_genome_fusion_jaffal/{sample}.txt'
    conda:
        '../envs/minimap2.yaml'
    shell:
        "minimap2 -x splice -t {threads} -c  {params.reference}  {input.fastq} -o {output.mapping_paf} 2>{log} "


rule minimap2_transcriptome_dup:
    input:
        transcriptome_file=config['reference']['transcriptome']
    output:
        transcriptome_dup='{project}/resources/transcriptome_dup.tsv'
    threads: config['threads']['minimap2']
    resources:
        mem_mb = 1024 * 30
    log:
        'logs/{project}/minimap2_transcriptome_dup.log'
    benchmark:
        'benchmarks/{project}/minimap2_transcriptome_dup.txt'
    conda:
        '../envs/minimap2.yaml'
    shell:
        "minimap2 {input.transcriptome_file} {input.transcriptome_file} -X -t {threads} -2 -c 2>{log} "
        "| cut -f1,6 "
        "| sed 's/_/\t/g' "
        "| awk 'BEGIN{{OFS=\"\\t\";}}{{print substr($1,1,15),substr($2,1,15),substr($3,1,15),substr($4,1,15);}}' "
        "| awk '$1!=$3' | sort | uniq > {output.transcriptome_dup} "





