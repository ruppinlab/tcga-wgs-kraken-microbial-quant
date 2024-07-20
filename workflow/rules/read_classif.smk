rule kraken2_db_taxonomy:
    params:
        db=KRAKEN2_DB_DIR,
        taskopt="--download-taxonomy",
        extra=config["kraken2"]["build"]["extra"],
    output:
        touch(KRAKEN2_DB_TAX_DONE_FILE),
    log:
        KRAKEN2_DB_TAX_LOG,
    threads: KRAKEN2_BUILD_THREADS
    wrapper:
        KRAKEN2_BUILD_WRAPPER


rule kraken2_db_library:
    params:
        db=KRAKEN2_DB_DIR,
        lib="{k2lib}",
        taskopt="--download-library",
        extra=config["kraken2"]["build"]["extra"],
    output:
        touch(KRAKEN2_DB_LIB_DONE_FILE),
    log:
        KRAKEN2_DB_LIB_LOG,
    threads: KRAKEN2_BUILD_THREADS
    wrapper:
        KRAKEN2_BUILD_WRAPPER


rule kraken2_db_build:
    input:
        KRAKEN2_DB_TAX_DONE_FILE,
        expand(KRAKEN2_DB_LIB_DONE_FILE, **EXPAND_PARAMS),
    params:
        db=KRAKEN2_DB_DIR,
        taskopt="--build",
        extra=config["kraken2"]["build"]["extra"],
    output:
        touch(KRAKEN2_DB_BUILD_DONE_FILE),
    log:
        KRAKEN2_DB_BUILD_LOG,
    threads: KRAKEN2_BUILD_THREADS
    wrapper:
        KRAKEN2_BUILD_WRAPPER


rule kraken2_read_classif:
    input:
        fqs=lambda wc: (
            [BOWTIE2_FILTERED_FASTQ_R1_FILE, BOWTIE2_FILTERED_FASTQ_R2_FILE]
            if wc.etype == "pe"
            else BOWTIE2_FILTERED_FASTQ_SE_FILE
        ),
        db=KRAKEN2_DB_DIR,
        build_done=KRAKEN2_DB_BUILD_DONE_FILE,
    params:
        extra=lambda wc: (
            f"{config['kraken2']['classify']['extra']['paired_end']} "
            f"{config['kraken2']['classify']['extra']['common']}"
            if wc.etype == "pe"
            else config["kraken2"]["classify"]["extra"]["common"]
        ),
    output:
        KRAKEN2_CLASSIF_FILE,
        report=KRAKEN2_REPORT_FILE,
    log:
        KRAKEN2_CLASSIF_LOG,
    threads: KRAKEN2_CLASSIFY_THREADS
    wrapper:
        KRAKEN2_CLASSIFY_WRAPPER


rule krakenuniq_db_taxonomy:
    params:
        db=KRAKENUNIQ_DB_DIR,
        lib="taxonomy",
        extra=config["krakenuniq"]["download"]["extra"],
    output:
        touch(KRAKENUNIQ_DB_TAX_DONE_FILE),
    log:
        KRAKENUNIQ_DB_TAX_LOG,
    threads: KRAKENUNIQ_DOWNLOAD_THREADS
    wrapper:
        KRAKENUNIQ_DOWNLOAD_WRAPPER


rule krakenuniq_db_library:
    params:
        db=KRAKENUNIQ_DB_DIR,
        lib="{kulib}",
        extra=config["krakenuniq"]["download"]["extra"],
    output:
        touch(KRAKENUNIQ_DB_LIB_DONE_FILE),
    log:
        KRAKENUNIQ_DB_LIB_LOG,
    threads: KRAKENUNIQ_DOWNLOAD_THREADS
    wrapper:
        KRAKENUNIQ_DOWNLOAD_WRAPPER


rule krakenuniq_db_build:
    input:
        KRAKENUNIQ_DB_TAX_DONE_FILE,
        expand(KRAKENUNIQ_DB_LIB_DONE_FILE, **EXPAND_PARAMS),
    params:
        db=KRAKENUNIQ_DB_DIR,
        taskopt="--build",
        extra=config["krakenuniq"]["build"]["extra"],
    output:
        touch(KRAKENUNIQ_DB_BUILD_DONE_FILE),
    log:
        KRAKENUNIQ_DB_BUILD_LOG,
    threads: KRAKENUNIQ_BUILD_THREADS
    wrapper:
        KRAKENUNIQ_BUILD_WRAPPER


rule krakenuniq_read_classif:
    input:
        fqs=lambda wc: (
            [BOWTIE2_FILTERED_FASTQ_R1_FILE, BOWTIE2_FILTERED_FASTQ_R2_FILE]
            if wc.etype == "pe"
            else BOWTIE2_FILTERED_FASTQ_SE_FILE
        ),
        db=KRAKENUNIQ_DB_DIR,
        build_done=KRAKENUNIQ_DB_BUILD_DONE_FILE,
    params:
        extra=lambda wc: (
            f"{config['krakenuniq']['classify']['extra']['paired_end']} "
            f"{config['krakenuniq']['classify']['extra']['common']}"
            if wc.etype == "pe"
            else config["krakenuniq"]["classify"]["extra"]["common"]
        ),
    output:
        KRAKENUNIQ_CLASSIF_FILE,
        report=KRAKENUNIQ_REPORT_FILE,
    log:
        KRAKENUNIQ_CLASSIF_LOG,
    threads: KRAKENUNIQ_CLASSIFY_THREADS
    wrapper:
        KRAKENUNIQ_CLASSIFY_WRAPPER
