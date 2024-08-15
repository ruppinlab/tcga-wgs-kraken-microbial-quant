rule kraken2_db_taxonomy:
    params:
        db=KRAKEN2_DB_DIR,
        task="download-taxonomy",
        protein=lambda wc: True if wc.k2dtype == "prot" else False,
        extra=config["kraken2"]["build"]["extra"],
    output:
        touch(KRAKEN2_DB_TAX_DONE_FILE),
    log:
        KRAKEN2_DB_TAX_LOG,
    wrapper:
        KRAKEN2_BUILD_WRAPPER


rule kraken2_nucl_db_library:
    input:
        KRAKEN2_NUCL_DB_TAX_DONE_FILE,
    params:
        k2=KRAKEN2_K2_SCRIPT_PATH,
        db=KRAKEN2_NUCL_DB_DIR,
        lib="{k2nlib}",
        task="download-library",
        protein=False,
        use_ftp=lambda wc: (
            True
            if config["kraken2"]["download"]["use_ftp"]
            or wc.k2nlib in config["resources"]["db"]["libs"]["kraken2"]["ftp_only"]
            else False
        ),
        extra=config["kraken2"]["k2"]["extra"],
        verbosity=config["kraken2"]["k2"]["verbosity"],
        backend=config["kraken2"]["k2"]["backend"],
        retries=config["download"]["retries"],
        retry_wait=config["download"]["retry_wait"],
    output:
        touch(KRAKEN2_NUCL_DB_LIB_DONE_FILE),
    log:
        KRAKEN2_NUCL_DB_LIB_LOG,
    threads: KRAKEN2_K2_THREADS
    wrapper:
        KRAKEN2_K2_WRAPPER


rule kraken2_prot_db_library:
    input:
        KRAKEN2_PROT_DB_TAX_DONE_FILE,
    params:
        k2=KRAKEN2_K2_SCRIPT_PATH,
        db=KRAKEN2_PROT_DB_DIR,
        lib="{k2plib}",
        task="download-library",
        protein=True,
        use_ftp=lambda wc: (
            True
            if config["kraken2"]["download"]["use_ftp"]
            or wc.k2plib in config["resources"]["db"]["libs"]["kraken2"]["ftp_only"]
            else False
        ),
        extra=config["kraken2"]["k2"]["extra"],
        verbosity=config["kraken2"]["k2"]["verbosity"],
        backend=config["kraken2"]["k2"]["backend"],
        retries=config["download"]["retries"],
        retry_wait=config["download"]["retry_wait"],
    output:
        touch(KRAKEN2_PROT_DB_LIB_DONE_FILE),
    log:
        KRAKEN2_PROT_DB_LIB_LOG,
    threads: KRAKEN2_K2_THREADS
    wrapper:
        KRAKEN2_K2_WRAPPER


rule kraken2_db:
    input:
        KRAKEN2_DB_TAX_DONE_FILE,
        lambda wc: expand(
            (
                KRAKEN2_NUCL_DB_LIB_DONE_FILE
                if wc.k2dtype == "nucl"
                else KRAKEN2_PROT_DB_LIB_DONE_FILE
            ),
            **EXPAND_PARAMS,
        ),
        lambda wc: (
            KRAKEN2_EUPATHDB_NUCL_LIB_IDMAP_FILE
            if wc.k2dtype == "nucl"
            else KRAKEN2_EUPATHDB_PROT_LIB_IDMAP_FILE
        ),
        lambda wc: (
            KRAKEN2_EUPATHDB_NUCL_LIB_FASTA_FILE
            if wc.k2dtype == "nucl"
            else KRAKEN2_EUPATHDB_PROT_LIB_FASTA_FILE
        ),
    params:
        db=KRAKEN2_DB_DIR,
        task="build",
        protein=lambda wc: True if wc.k2dtype == "prot" else False,
        extra=config["kraken2"]["build"]["extra"],
    output:
        touch(KRAKEN2_DB_DONE_FILE),
    log:
        KRAKEN2_DB_LOG,
    threads: KRAKEN2_BUILD_THREADS
    wrapper:
        KRAKEN2_BUILD_WRAPPER


rule krakenuniq_db_taxonomy:
    params:
        db=KRAKENUNIQ_DB_DIR,
        lib="taxonomy",
        rsync=True,
        extra=config["krakenuniq"]["download"]["extra"],
    output:
        touch(KRAKENUNIQ_DB_TAX_DONE_FILE),
    log:
        KRAKENUNIQ_DB_TAX_LOG,
    threads: KRAKENUNIQ_BUILD_THREADS
    wrapper:
        KRAKENUNIQ_DOWNLOAD_WRAPPER


rule krakenuniq_db_library:
    input:
        KRAKENUNIQ_DB_TAX_DONE_FILE,
    params:
        db=KRAKENUNIQ_DB_DIR,
        lib="{kulib}",
        rsync=True if config["krakenuniq"]["download"]["use_rsync"] else False,
        extra=config["krakenuniq"]["download"]["extra"],
    output:
        touch(KRAKENUNIQ_DB_LIB_DONE_FILE),
    log:
        KRAKENUNIQ_DB_LIB_LOG,
    threads: KRAKENUNIQ_BUILD_THREADS
    wrapper:
        KRAKENUNIQ_DOWNLOAD_WRAPPER


rule krakenuniq_db:
    input:
        KRAKENUNIQ_DB_TAX_DONE_FILE,
        expand(KRAKENUNIQ_DB_LIB_DONE_FILE, **EXPAND_PARAMS),
    params:
        db=KRAKENUNIQ_DB_DIR,
        task="build",
        extra=config["krakenuniq"]["build"]["extra"],
    output:
        touch(KRAKENUNIQ_DB_DONE_FILE),
    log:
        KRAKENUNIQ_DB_LOG,
    threads: KRAKENUNIQ_BUILD_THREADS
    wrapper:
        KRAKENUNIQ_BUILD_WRAPPER


rule bracken_db_kraken_classif:
    input:
        db_done=(
            KRAKEN2_NUCL_DB_DONE_FILE
            if KRAKEN_MODE == "kraken2"
            else KRAKENUNIQ_DB_DONE_FILE
        ),
    params:
        bracken_build=BRACKEN_BUILD_SCRIPT_PATH,
        db=KRAKEN2_NUCL_DB_DIR if KRAKEN_MODE == "kraken2" else KRAKENUNIQ_DB_DIR,
        klen=35 if KRAKEN_MODE == "kraken2" else 31,
        ktype=KRAKEN_MODE,
        readlen="{readlen}",
        db_only=True,
    output:
        touch(BRACKEN_DB_KRAKEN_CLASSIF_DONE_FILE),
    log:
        BRACKEN_DB_KRAKEN_CLASSIF_LOG,
    threads: BRACKEN_BUILD_THREADS
    wrapper:
        BRACKEN_BUILD_WRAPPER


rule bracken_db_kmer_distr:
    input:
        db_done=BRACKEN_DB_KRAKEN_CLASSIF_DONE_FILE,
    params:
        bracken_build=BRACKEN_BUILD_SCRIPT_PATH,
        db=KRAKEN2_NUCL_DB_DIR if KRAKEN_MODE == "kraken2" else KRAKENUNIQ_DB_DIR,
        klen=35 if KRAKEN_MODE == "kraken2" else 31,
        ktype=KRAKEN_MODE,
        readlen="{readlen}",
    output:
        touch(BRACKEN_DB_KMER_DISTR_DONE_FILE),
    log:
        BRACKEN_DB_KMER_DISTR_LOG,
    threads: BRACKEN_BUILD_THREADS
    wrapper:
        BRACKEN_BUILD_WRAPPER
