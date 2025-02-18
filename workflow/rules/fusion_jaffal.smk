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
        "process_transcriptome_align_table {input.mapping_paf} 1000  {transtable} > {output.gene_count} 2>{log}"

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
    shell:
        "cat {input.gene_count} | awk '{print \$1}' > {output.fusion_temp} && "
        "reformat.sh in={input.raw_fasta} out=stdout.fasta fastawrap=0 "
        "| extract_seq_from_fasta {output.fusion_temp} > {output.fusion_fa}"


# minimap2_genome = {
#    doc "Aligning candidates to genome using minimap2"
#    output.dir=jaffa_output+branch
#    produce(branch+"_genome.paf",branch+"_genome.psl"){ 
# 	exec """
# 	   $minimap2 -t $threads -x splice -c $genomeFasta $input > $output1;
# 	   grep \$'\\t+\\t' $output1 | awk -F'\\t' -v OFS="\\t" '{ print \$4-\$3,0,0,0,0,0,0,0,\$5,\$1,\$2,\$3,\$4,\$6,\$7,\$8,\$9,2, 100","\$4-\$3-100",",\$3","\$3+100",",  \$8","\$9-\$4+\$3+100"," }' > $output2 ;
# 	   grep \$'\\t-\\t' $output1 | awk -F'\\t' -v OFS="\\t" '{ print \$4-\$3,0,0,0,0,0,0,0,\$5,\$1,\$2,\$3,\$4,\$6,\$7,\$8,\$9,2, 100","\$4-\$3-100",", \$2-\$4","\$2-\$4+100",", \$8","\$9-\$4+\$3+100"," }' >> $output2 ;
#         """
#    }
# }

# bpipe="/opt/bin/bpipe"
# samtools="/opt/bin/samtools"
# R="/usr/bin/R"
# minimap2="/opt/bin/minimap2"
# process_transcriptome_align_table="/opt/bin/process_transcriptome_align_table"
# make_simple_read_table="/opt/bin/make_simple_read_table"
# extract_seq_from_fasta="/opt/bin/extract_seq_from_fasta"
# reformat="/opt/bin/reformat.sh"
# make_3_gene_fusion_table="/opt/bin/make_3_gene_fusion_table"


# rule jaffal_fusion:
#     input:
#         mapping_bam='{project}/{sample}/alignment/{sample}_minimap2_genome_4longgf.bam',
#     output:
#         fusions = '{project}/{sample}/fusion/longgf/{sample}_fusion_longgf.tsv'
#     params:
#         gtf=config['reference']['annotation']
#     log:
#         'logs/{project}/jaffal_fusion/{sample}.log'
#     benchmark:
#         'benchmarks/{project}/jaffal_fusion/{sample}.txt'
#     resources:
#         mem_mb = 1024 * 10
#     threads: 1
#     container:
#         "docker://btrspg/longgf:0.1.2"
#     shell:
#         'LongGF {input.mapping_bam}  {params.gtf} 100 30 100 '
#         '| grep "SumGF" > {output.fusions}'


