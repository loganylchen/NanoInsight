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


rule vg_genome_index_aeron:
    input:
        reference = config['reference']['genome'],
        annotation= config['reference']['annotation'],
    output:
        '{project}/resources/genome.gfa',
    log:
        'logs/{project}/vg_genome_index_aeron.log'
    container:
        'docker://btrspg/aeron:a6e7d589e3feeb22b5374b455a1a677e3bb2edfa'
    threads: 1
    benchmark:
        'benchmarks/{project}/vg_genome_index_aeron.txt'
    resources:
        mem_mb = 1024 * 300,
    shell:
        "GraphBuilder.py -e {input.reference} -g {input.annotation} -o {output} 2>{log} 1>&2  "