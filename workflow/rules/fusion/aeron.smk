rule partial_pairs_aeron:
    input:
        readalin = '{project}/{sample}/alignment/{sample}_vg_genome_4aeron.gam',
        reads = get_raw_fastq
    output: 
        temp('{project}/{sample}/alignment/{sample}_pairs_4aeron.gam')
    log:
        'logs/{project}/fusion_preprocessing_pair_aeron/{sample}.log'
    benchmark:
        'benchmarks/{project}/fusion_preprocessing_pair_aeron/{sample}.txt' 
    container:
        'docker://btrspg/aeron:a6e7d589e3feeb22b5374b455a1a677e3bb2edfa'
    threads: config['threads']['aeron']
    resources:
        mem_mb = 1024 * 30,
    shell: 
        "PickAdjacentAlnPairs {input.readalin} 20 {input.reads} {output} 100 2>{log}"


rule postprocess_aeron:
    input:
        all_alns = '{project}/{sample}/alignment/{sample}_graphaligner_genome_4aeron.gam',
        reads = get_raw_fastq
    output:
        selected_alns = '{project}/{sample}/alignment/{sample}_selected_4aeron.gam',
        full_len_alns = '{project}/{sample}/alignment/{sample}_fulllength_4aeron.gam',
        summary = '{project}/{sample}/alignment/{sample}_4aeron.summary.txt',
    container:
        'docker://btrspg/aeron:a6e7d589e3feeb22b5374b455a1a677e3bb2edfa'
    log:
        'logs/{project}/postprocess_fusion_aeron/{sample}.log'
    benchmark:
        'benchmarks/{project}/postprocess_fusion_aeron/{sample}.txt'
    threads: config['threads']['aeron']
    shell:
        "Postprocess {input.all_alns} {input.reads} {output.selected_alns} {output.full_len_alns} {output.summary} 2>{log} 1>&2"

rule pair_assignments_aeron:
    input:
        transcriptaln = '{project}/{sample}/alignment/{sample}_fulllength_4aeron.gam',
        pairs = '{project}/{sample}/alignment/{sample}_pairs_4aeron.gam',
        reads = get_raw_fastq
    output: 
        temp('{project}/{sample}/fusion/aeron/{sample}_pairs_4aeron.txt')
    container:
        'docker://btrspg/aeron:a6e7d589e3feeb22b5374b455a1a677e3bb2edfa'
    log:
        'logs/{project}/fusion_preprocessing_pair_aeron/{sample}.log'
    benchmark:
        'benchmarks/{project}/fusion_preprocessing_pair_aeron/{sample}.txt'
    threads: config['threads']['aeron']
    shell: 
        "AlignmentSubsequenceIdentity {input.transcriptaln} {input.pairs} {input.reads} 1 > {output}"



rule pick_longest_aeron:
    input:
        '{project}/{sample}/alignment/{sample}_selected_4aeron.gam',
    output:
        '{project}/{sample}/alignment/{sample}_longest_4aeron.gam',
    benchmark:
        "benchmarks/{project}/pick_longest_{sample}.txt"
    log:
        "logs/{project}/pick_longest_aeron/{sample}.log"
    threads: config['threads']['aeron']
    container:
        'docker://btrspg/aeron:a6e7d589e3feeb22b5374b455a1a677e3bb2edfa'
    shell:
        "SelectLongestAlignment {input} {output} 2>{log} 1>&2"


rule exact_pair_assignments_aeron:
    input: 
        "{project}/{sample}/fusion/aeron/{sample}_pairs_4aeron.txt"
    output: 
        temp("{project}/{sample}/fusion/aeron/{sample}_exactMatrix.txt")
    benchmark:
        "benchmarks/{project}/exact_pair_assignments_{sample}.txt"
    log:
        "logs/{project}/exact_pair_assignments_aeron/{sample}.log"
    threads: config['threads']['aeron']
    container:
        'docker://btrspg/aeron:a6e7d589e3feeb22b5374b455a1a677e3bb2edfa'
    shell: 
        "awk -F '\\t' '{{if ($3 == 1) print;}}' < {input} > {output}"

rule parse_matrix_aeron:
    input: 
        "{project}/{sample}/fusion/aeron/{sample}_exactMatrix.txt"
    output: 
        temp("{project}/{sample}/fusion/aeron/{sample}_exactMatrixparse.txt")
    benchmark:
        "benchmarks/{project}/parse_exact_pair_assignments_{sample}.txt"
    log:
        "logs/{project}/parse_exact_pair_assignments_aeron/{sample}.log"
    threads: config['threads']['aeron']
    container:
        'docker://btrspg/aeron:a6e7d589e3feeb22b5374b455a1a677e3bb2edfa'
    shell: 
        "parse_pair_matrix.py < {input} > {output}"

rule loose_fusions_aeron:
    input: 
        "{project}/{sample}/fusion/aeron/{sample}_exactMatrixparse.txt"
    output: 
        temp("{project}/{sample}/fusion/aeron/{sample}_loose_fusions.txt")
    benchmark:
        "benchmarks/{project}/loose_fusions_{sample}.txt"
    log:
        "logs/{project}/loose_fusions_aeron/{sample}.log"
    threads: config['threads']['aeron']
    container:
        'docker://btrspg/aeron:a6e7d589e3feeb22b5374b455a1a677e3bb2edfa'
    shell: 
        "pairmatrix_get_genes.py < {input} > {output}"



rule fusionfinder_aeron:
    input:
        graph = '{project}/resources/genome.gfa',
        loose_fusions = "{project}/{sample}/fusion/aeron/{sample}_loose_fusions.txt",
        pairmatrix = "{project}/{sample}/fusion/aeron/{sample}_exactMatrix.txt",
        transcriptaln = '{project}/{sample}/alignment/{sample}_fulllength_4aeron.gam',
        reads = get_raw_fastq
    output:
        fusions = "{project}/{sample}/fusion/aeron/{sample}_unfiltered_fusion.txt",
        corrected ="{project}/{sample}/fusion/aeron/{sample}_corrected_fusion.txt",
    log:
        "logs/{project}/fusionfinder_aeron/{sample}.log"
    benchmark:
        "benchmarks/{project}/fusionfinder_aeron/{sample}.txt"
    threads: config['threads']['aeron']
    container:
        'docker://btrspg/aeron:a6e7d589e3feeb22b5374b455a1a677e3bb2edfa'
    shell: 
        "FusionFinder {input.graph} {input.loose_fusions} {input.pairmatrix} {input.transcriptaln} {input.reads} "
        " 1 1.0 1 1 {threads} {output.fusions} {output.corrected} 2>{log} 1>&2"

rule filter_fusions_aeron:
    input: 
        "{project}/{sample}/fusion/aeron/{sample}_unfiltered_fusion.txt",
    output: 
        "{project}/{sample}/fusion/aeron/{sample}_fusion_table.txt",
    log:
        "logs/{project}/filter_fusions_aeron/{sample}.log"
    benchmark:
        "benchmarks/{project}/filter_fusions_aeron/{sample}.txt"
    threads: config['threads']['aeron']
    container:
        'docker://btrspg/aeron:a6e7d589e3feeb22b5374b455a1a677e3bb2edfa'
    params:
        FUSION_MAX_ERROR_RATE = 0.2,
        FUSION_MIN_SCORE_DIFFERENCE = 200
    shell: 
        "awk -F '\t' '{{if ($2 < {params.FUSION_MAX_ERROR_RATE} && $3 < -{params.FUSION_MIN_SCORE_DIFFERENCE}) print;}}' < {input} > {output}"

rule fusion_transcripts_aeron:
    input:
        fusions = "{project}/{sample}/fusion/aeron/{sample}_fusion_table.txt",
        corrected = "{project}/{sample}/fusion/aeron/{sample}_corrected_fusion.txt",
    output: 
        "{project}/{sample}/fusion/aeron/{sample}.fa",
    log:
        "logs/{project}/fusion_transcripts_aeron/{sample}.log"
    benchmark:
        "benchmarks/{project}/fusion_transcripts_aeron/{sample}.txt"
    threads: config['threads']['aeron']
    container:
        'docker://btrspg/aeron:a6e7d589e3feeb22b5374b455a1a677e3bb2edfa'
    shell: "pick_fusion_exemplar.py {input.fusions} {input.corrected} > {output} 2>{log}"

rule merge_ref_and_fusions_aeron:
    input:
        "{project}/{sample}/fusion/aeron/{sample}.fa",
    params:
        reference = config['reference']['transcriptome']
    output: 
        "{project}/{sample}/fusion/aeron/{sample}_withref.fa",
    threads: 1
    shell: 
        "cat {params.reference} {input.fusions} > {output}"

rule align_reads_to_fusions_aeron:
    input:
        refplusfusion = "{project}/{sample}/fusion/aeron/{sample}_withref.fa",
        reads = get_raw_fastq
    output: 
        temp("{project}/{sample}/fusion/aeron/{sample}_tofusion.bam")
    log:
        "logs/{project}/align_reads_to_fusions_aeron/{sample}.log"    
    threads: config['threads']['minimap2']
    container:
        'docker://btrspg/minimap2:2.28'
    benchmark:
        "benchmarks/{project}/align_reads_to_fusions_aeron/{sample}.txt"
    shell: "minimap2 --secondary=no -t {threads} -a -x map-ont {input.refplusfusion} {input.reads} "
    "| samtools view -bSh "
    "| samtools sort - -o {output} 2>{log} && "
    "samtools index {output}"



rule fusion_support_sam_aeron:
    input:
        fusiontranscripts = "{project}/{sample}/fusion/aeron/{sample}.fa",
        minimapalns = "{project}/{sample}/fusion/aeron/{sample}_tofusion.bam"
    output:
        supportfile = "{project}/{sample}/fusion/aeron/{sample}_fusion_support.txt",
        fusionalnfile = "{project}/{sample}/fusion/aeron/{sample}_tofusion.sam"
    log:
        "logs/{project}/fusion_support_sam_aeron/{sample}.log"
    benchmark:
        "benchmarks/{project}/fusion_support_sam_aeron/{sample}.txt"
    threads: config['threads']['aeron']
    container:
        'docker://btrspg/aeron:a6e7d589e3feeb22b5374b455a1a677e3bb2edfa'
    shell: "fusion_support_sam.py {input.fusiontranscripts} {input.minimapalns} 150 {output.supportfile} {output.fusionalnfile}"











