rule eupathdb_fasta_archive:
    params:
        EUPATHDB_FASTA_TGZ_URL,
    output:
        temp(EUPATHDB_FASTA_TGZ_FILE),
    log:
        EUPATHDB_FASTA_TGZ_LOG,
    message:
        "{params}"
    retries: config["download"]["retries"]
    conda:
        "../envs/wget.yaml"
    shell:
        "wget -nv -O {output} {params} > {log} 2>&1"


rule eupathdb_fastas:
    input:
        EUPATHDB_FASTA_TGZ_FILE,
    params:
        EUPATHDB_RESOURCES_DIR,
    output:
        temp(directory(EUPATHDB_FASTA_DIR)),
    log:
        EUPATHDB_FASTA_LOG,
    shell:
        "tar -xvzf {input} -C {params} > {log} 2>&1"


rule eupathdb_merged_fasta:
    input:
        EUPATHDB_FASTA_DIR,
    output:
        EUPATHDB_MERGED_FASTA_FILE,
    log:
        EUPATHDB_MERGED_FASTA_LOG,
    shell:
        "find {input} -type f -name '*.fna' -exec cat {{}} + 1> {output} 2> {log}"


rule eupathdb_seqid2taxid_map:
    params:
        EUPATHDB_SEQID2TAXID_MAP_URL,
    output:
        temp(EUPATHDB_SEQID2TAXID_MAP_FILE),
    log:
        EUPATHDB_SEQID2TAXID_MAP_LOG,
    message:
        "{params}"
    retries: config["download"]["retries"]
    conda:
        "../envs/wget.yaml"
    shell:
        "wget -nv -O {output} {params} > {log} 2>&1"


rule kraken2_eupathdb_library:
    input:
        idmap=EUPATHDB_SEQID2TAXID_MAP_FILE,
        fasta=EUPATHDB_MERGED_FASTA_FILE,
    output:
        idmap=KRAKEN2_EUPATHDB_LIB_IDMAP_FILE,
        fasta=KRAKEN2_EUPATHDB_LIB_FASTA_FILE,
    log:
        KRAKEN2_EUPATHDB_LIB_FASTA_LOG,
    script:
        "../scripts/kraken2_eupathdb_library.py"
