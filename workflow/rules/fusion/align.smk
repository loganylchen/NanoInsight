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
    container:
        'docker://btrspg/minimap2:2.28'
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
    priority: 1
    resources:
        mem_mb = 1024 * 30
    log:
        'logs/{project}/minimap2_genome_fusion_longgf/{sample}.log'
    benchmark:
        'benchmarks/{project}/minimap2_genome_fusion_longgf/{sample}.txt'
    container:
        'docker://btrspg/minimap2:2.28'
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
    priority: 1
    resources:
        mem_mb = 1024 * 30
    log:
        'logs/{project}/minimap2_transcriptome_fusion_jaffal/{sample}.log'
    benchmark:
        'benchmarks/{project}/minimap2_transcriptome_fusion_jaffal/{sample}.txt'
    container:
        'docker://btrspg/minimap2:2.28'
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
    priority: 1
    resources:
        mem_mb = 1024 * 30
    log:
        'logs/{project}/minimap2_genome_fusion_jaffal/{sample}.log'
    benchmark:
        'benchmarks/{project}/minimap2_genome_fusion_jaffal/{sample}.txt'
    container:
        'docker://btrspg/minimap2:2.28'
    shell:
        "minimap2 -x splice -t {threads} -c  {params.reference}  {input.fastq} -o {output.mapping_paf} 2>{log} "


rule minimap2_transcriptome_dup:
    input:
        transcriptome_file=config['reference']['transcriptome']
    output:
        transcriptome_dup='{project}/resources/transcriptome_dup.tsv'
    threads: config['threads']['minimap2']
    priority: 1
    resources:
        mem_mb = 1024 * 30
    log:
        'logs/{project}/minimap2_transcriptome_dup.log'
    benchmark:
        'benchmarks/{project}/minimap2_transcriptome_dup.txt'
    container:
        'docker://btrspg/minimap2:2.28'
    shell:
        "minimap2 {input.transcriptome_file} {input.transcriptome_file} -X -t {threads} -2 -c 2>{log} "
        "| cut -f1,6 "
        "| sed 's/_/\t/g' "
        "| awk 'BEGIN{{OFS=\"\\t\";}}{{print substr($1,1,15),substr($2,1,15),substr($3,1,15),substr($4,1,15);}}' "
        "| awk '$1!=$3' | sort | uniq > {output.transcriptome_dup} "


rule vg_gfa:
    input:
        reference=config['reference']['genome'],
    output:
        vg_gfa='{project}/resources/genome.gfa',
        vg_vg = '{project}/resources/genome.vg',
    threads: config['threads']['vg']
    log:
        'logs/{project}/vg_gfa.log'
    benchmark:
        'benchmarks/{project}/vg_gfa.txt'
    container:
        'docker://btrspg/vg:1.23.0'
    shell:
        "vg construct -r {input.reference} -t {threads} > {output.vg_vg} 2>{log} && "
        "vg view {output.vg_vg} > {output.vg_gfa} 2>>{log}"

rule vg_genome_aeron:
    input:
        graph = '{project}/resources/genome.gfa',
        reads = get_raw_fastq
    output:
        '{project}/{sample}/alignment/{sample}_vg_genome_4aeron.gam'
    log:
        'logs/{project}/vg_genome_fusion_aeron/{sample}.log'
    params:
        seedsize = 17,
        maxseeds=15,
        aligner_bandwidth = 35
    container:
        'docker://btrspg/aeron:a6e7d589e3feeb22b5374b455a1a677e3bb2edfa'
    threads: config['threads']['vg']
    benchmark:
        'benchmarks/{project}/vg_genome_fusion_aeron/{sample}.txt'
    resources:
        mem_mb = 1024 * 30,
        tmpdir = '{project}/{sample}/alignment/{sample}_vg_tmp'
    shell:
        "Aligner --all-alignments -g {input.graph} -f {input.reads} "
        "--try-all-seeds --seeds-mxm-length {params.seedsize} --seeds-mem-count {params.maxseeds} "
        "--seeds-mxm-cache-prefix {resources.tmpdir} "
        "-a {output} -t {threads} -b {params.aligner_bandwidth} --E-cutoff 1  2> {log} 1>&2 "

rule align_aeron:
    input:
        graph = '{project}/resources/genome.gfa',
        reads = get_raw_fastq
    output:
        '{project}/{sample}/alignment/{sample}_graphaligner_genome_4aeron.gam'
    benchmark:
        'benchmarks/{project}/graphaligner_fusion_aeron/{sample}.txt'
    log:
        'logs/{project}/graphaligner_fusion_aeron/{sample}.log'
    threads: config['threads']['aeron']
    container:
        'docker://btrspg/aeron:a6e7d589e3feeb22b5374b455a1a677e3bb2edfa'
    params:
        seedsize = 17,
        maxseeds=15,
        aligner_bandwidth = 35
    resources:
        mem_mb = 1024 * 30,
        tmpdir = '{project}/{sample}/alignment/{sample}_ga_tmp'
    shell:
        "GraphAligner -g {input.graph} "
        "-f {input.reads} --try-all-seeds --seeds-mxm-length {params.seedsize} 
        "--seeds-mem-count {params.maxseeds} "
        "--seeds-mxm-cache-prefix {resources.tmpdir} "
        "-a {output} -t {threads} -b {params.aligner_bandwidth} --greedy-length 2>{log} 1>&2"

