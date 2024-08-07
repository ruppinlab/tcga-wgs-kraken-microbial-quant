# Bioconductor EuPathDB library doesn't produce as comprehesnive species
# and taxonomy metadata as manual downloading from the web, metadata is
# different and misses a lot of organisms
# rule eupathdb_metadata:
#     output:
#         EUPATHDB_METADATA_FILE,
#     log:
#         EUPATHDB_METADATA_LOG,
#     conda:
#         "../envs/eupathdb.yaml"
#     script:
#         "../scripts/eupathdb_metadata.R"


rule eupathdb_fasta_archive:
    params:
        lambda wc: (
            EUPATHDB_NUCL_FASTA_TGZ_URL
            if wc.k2dtype == "nucl"
            else EUPATHDB_PROT_FASTA_TGZ_URL
        ),
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
    run:
        shell("tar -xvzf {input} -C {params} > {log} 2>&1")
        if wildcards.k2dtype == "nucl":
            shell(f"mv -vf {EUPATHDB_NUCL_RESOURCES_DIR}/eupathDB54_CLEAN {output}")


rule eupathdb_merged_fasta:
    input:
        EUPATHDB_FASTA_DIR,
    output:
        temp(EUPATHDB_MERGED_FASTA_FILE),
    log:
        EUPATHDB_MERGED_FASTA_LOG,
    shell:
        "find {input} -type f "
        "-regextype awk -regex '.+?\\.(fna|fasta|faa|fa)$' "
        "-exec cat {{}} + 1> {output} 2> {log}"


rule eupathdb_nucl_seqid2taxid_map:
    params:
        EUPATHDB_NUCL_SEQID2TAXID_MAP_URL,
    output:
        temp(EUPATHDB_NUCL_SEQID2TAXID_MAP_FILE),
    log:
        EUPATHDB_NUCL_SEQID2TAXID_MAP_LOG,
    message:
        "{params}"
    retries: config["download"]["retries"]
    conda:
        "../envs/wget.yaml"
    shell:
        "wget -nv -O {output} {params} > {log} 2>&1"


rule kraken2_eupathdb_nucl_library:
    input:
        idmap=EUPATHDB_NUCL_SEQID2TAXID_MAP_FILE,
        fasta=EUPATHDB_NUCL_MERGED_FASTA_FILE,
    output:
        idmap=KRAKEN2_EUPATHDB_NUCL_LIB_IDMAP_FILE,
        fasta=KRAKEN2_EUPATHDB_NUCL_LIB_FASTA_FILE,
    log:
        KRAKEN2_EUPATHDB_NUCL_LIB_FASTA_LOG,
    script:
        "../scripts/kraken2_eupathdb_nucl_library.py"


rule kraken2_eupathdb_prot_library:
    input:
        meta=EUPATHDB_METADATA_FILE,
        fasta=EUPATHDB_PROT_MERGED_FASTA_FILE,
    output:
        idmap=KRAKEN2_EUPATHDB_PROT_LIB_IDMAP_FILE,
        fasta=KRAKEN2_EUPATHDB_PROT_LIB_FASTA_FILE,
    log:
        KRAKEN2_EUPATHDB_PROT_LIB_FASTA_LOG,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/kraken2_eupathdb_prot_library.py"
