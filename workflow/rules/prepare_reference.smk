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



