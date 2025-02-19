rule jaffal_fastq_to_fasta:
    input:
        raw_fastq=get_raw_fastq    
    output:
        fasta='{project}/{sample}/fusion/jaffal/raw.fasta'
    threads: 1
    resources:
        mem_mb = 1024 * 5
    benchmark:
        'benchmarks/{project}/jaffal_fastq_to_fasta/{sample}.txt'
    log:
        'logs/{project}/jaffal_fastq_to_fasta/{sample}.log'
    shell:
        "seqtk seq -a {input.raw_fastq} > {output.fasta} 2>{log}"

rule jaffal_filter_genome_mapping:
    input:
        mapping_paf='{project}/{sample}/alignment/{sample}_minimap2_genome_4jaffal.paf'
    output:
        genome_psl='{project}/{sample}/fusion/jaffal/{sample}.genome.psl'
    log:
        'logs/{project}/jaffal_filter_genome_mapping/{sample}.log'
    benchmark:
        'benchmarks/{project}/jaffal_filter_genome_mapping/{sample}.txt'
    resources:
        mem_mb= 1024*10
    threads: 1
    script:
        '../scripts/jaffal_generate_genome_psl.sh'



rule jaffal_filter_transcripts:
    input:
        mapping_paf='{project}/{sample}/alignment/{sample}_minimap2_transcriptome_4jaffal.paf'
    output:
        gene_count='{project}/{sample}/fusion/jaffal/{sample}.geneCounts.tsv'
    params:
        transtable=config['reference']['transtable']
    log:
        'logs/{project}/jaffal_filter_transcripts/{sample}.log'
    benchmark:
        'benchmarks/{project}/jaffal_filter_transcripts/{sample}.txt'
    resources:
        mem_mb= 1024*10
    threads: 1
    container:
        "docker://btrspg/jaffal:2.3"
    shell:
        "process_transcriptome_align_table {input.mapping_paf} 1000  {params.transtable} > {output.gene_count} 2>{log}"

rule jaffal_extract_sequences:
    input:
        gene_count='{project}/{sample}/fusion/jaffal/{sample}.geneCounts.tsv',
        raw_fasta='{project}/{sample}/fusion/jaffal/raw.fasta'
    output:
        fusion_fa='{project}/{sample}/fusion/jaffal/{sample}.fusion.fa',
        fusion_temp=temp('{project}/{sample}/fusion/jaffal/{sample}.temp')
    log:
        'logs/{project}/jaffal_extract_sequences/{sample}.log'
    benchmark:
        'benchmarks/{project}/jaffal_extract_sequences/{sample}.txt'
    resources:
        mem_mb= 1024*10
    threads: 1
    container:
        "docker://btrspg/jaffal:2.3"
    script:
        '../scripts/jaffal_extract_sequences.sh'


rule jaffal_make_fasta_reads_table:
    input:
        gene_count='{project}/{sample}/fusion/jaffal/{sample}.geneCounts.tsv', 
    output:
        fasta_reads_table='{project}/{sample}/fusion/jaffal/{sample}.geneCounts.reads'
    log:
        'logs/{project}/jaffal_make_fasta_reads_table/{sample}.log'
    benchmark:
        'benchmarks/{project}/jaffal_make_fasta_reads_table/{sample}.txt'
    resources:
        mem_mb= 1024*10
    threads: 1
    container:
        "docker://btrspg/jaffal:2.3"
    script:
        '../scripts/jaffal_make_fasta_reads_table.sh'



rule jaffal_get_final_list:
    input:
        genome_psl='{project}/{sample}/fusion/jaffal/{sample}.genome.psl',
        fasta_reads_table='{project}/{sample}/fusion/jaffal/{sample}.geneCounts.reads',
        gene_count='{project}/{sample}/fusion/jaffal/{sample}.geneCounts.tsv', 
    output:
        final_list='{project}/{sample}/fusion/jaffal/{sample}.final_list.txt',
        three_gene_reads='{project}/{sample}/fusion/jaffal/{sample}.3gene.reads',
        three_gene_summary='{project}/{sample}/fusion/jaffal/{sample}.3gene.summary.txt',
    log:
        'logs/{project}/jaffal_get_final_list/{sample}.log'
    params:
        transtable=config['reference']['transtable'],
        knowntable=config['reference']['knowntable'],
    benchmark:
        'benchmarks/{project}/jaffal_get_final_list/{sample}.txt'
    resources:
        mem_mb= 1024*10
    threads: 1
    container:
        "docker://btrspg/jaffal:2.3"
    script:
        '../scripts/jaffal_get_final_list.sh'

rule jaffal_compile_results:
    input:
        three_gene_summary='{project}/{sample}/fusion/jaffal/{sample}.3gene.summary.txt',
        final_list='{project}/{sample}/fusion/jaffal/{sample}.final_list.txt',
         
    output:
        output_fasta='{project}/{sample}/fusion/jaffal/{sample}.final.fasta',
        fusion='{project}/{sample}/fusion/jaffal/{sample}_fusion_jaffal.tsv',
    log:
        'logs/{project}/jaffal_compile_results/{sample}.log'
    params:
        final_out_prefix='{project}/{sample}/fusion/jaffal/{sample}.final'
    benchmark:
        'benchmarks/{project}/jaffal_compile_results/{sample}.txt'
    resources:
        mem_mb= 1024*10
    threads: 1
    container:
        "docker://btrspg/jaffal:2.3"
    script:
        '../scripts/jaffal_compile_results.sh'