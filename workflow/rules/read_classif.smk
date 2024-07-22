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
    threads: KRAKEN2_BUILD_THREADS
    wrapper:
        KRAKEN2_BUILD_WRAPPER


rule kraken2_nucl_library:
    params:
        db=KRAKEN2_NUCL_DB_DIR,
        lib="{k2nlib}",
        task="download-library",
        protein=False,
        extra=config["kraken2"]["build"]["extra"],
    output:
        touch(KRAKEN2_NUCL_DB_LIB_DONE_FILE),
    log:
        KRAKEN2_NUCL_DB_LIB_LOG,
    threads: KRAKEN2_BUILD_THREADS
    wrapper:
        KRAKEN2_BUILD_WRAPPER


rule kraken2_prot_library:
    params:
        db=KRAKEN2_PROT_DB_DIR,
        lib="{k2plib}",
        task="download-library",
        protein=True,
        extra=config["kraken2"]["build"]["extra"],
    output:
        touch(KRAKEN2_PROT_DB_LIB_DONE_FILE),
    log:
        KRAKEN2_PROT_DB_LIB_LOG,
    threads: KRAKEN2_BUILD_THREADS
    wrapper:
        KRAKEN2_BUILD_WRAPPER


rule kraken2_db_build:
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
    params:
        db=KRAKEN2_DB_DIR,
        task="build",
        protein=lambda wc: True if wc.k2dtype == "prot" else False,
        extra=config["kraken2"]["build"]["extra"],
    output:
        touch(KRAKEN2_DB_BUILD_DONE_FILE),
    log:
        KRAKEN2_DB_BUILD_LOG,
    threads: KRAKEN2_BUILD_THREADS
    wrapper:
        KRAKEN2_BUILD_WRAPPER


rule kraken2_nucl_read_classif_pe:
    input:
        fqs=[BOWTIE2_FILTERED_FASTQ_R1_FILE, BOWTIE2_FILTERED_FASTQ_R2_FILE],
        db=KRAKEN2_NUCL_DB_DIR,
        build_done=KRAKEN2_NUCL_DB_BUILD_DONE_FILE,
    params:
        output="-",
        paired_end=True,
        extra=(
            f"{config['kraken2']['classify']['extra']['paired_end']} "
            f"{config['kraken2']['classify']['extra']['common']}"
        ),
    output:
        classif=[
            temp(KRAKEN2_NUCL_CLASSIF_FASTQ_R1_DC_FILE),
            temp(KRAKEN2_NUCL_CLASSIF_FASTQ_R2_DC_FILE),
        ],
        unclassif=[
            temp(KRAKEN2_NUCL_UNCLASSIF_FASTQ_R1_DC_FILE),
            temp(KRAKEN2_NUCL_UNCLASSIF_FASTQ_R2_DC_FILE),
        ],
        report=KRAKEN2_NUCL_REPORT_PE_FILE,
    log:
        KRAKEN2_NUCL_CLASSIFY_PE_LOG,
    threads: KRAKEN2_CLASSIFY_THREADS
    wrapper:
        KRAKEN2_CLASSIFY_WRAPPER


rule kraken2_nucl_read_classif_se:
    input:
        fqs=BOWTIE2_FILTERED_FASTQ_SE_FILE,
        db=KRAKEN2_NUCL_DB_DIR,
        build_done=KRAKEN2_NUCL_DB_BUILD_DONE_FILE,
    params:
        output="-",
        paired_end=False,
        extra=config["kraken2"]["classify"]["extra"]["common"],
    output:
        classif=temp(KRAKEN2_NUCL_CLASSIF_FASTQ_SE_DC_FILE),
        unclassif=temp(KRAKEN2_NUCL_UNCLASSIF_FASTQ_SE_DC_FILE),
        report=KRAKEN2_NUCL_REPORT_SE_FILE,
    log:
        KRAKEN2_NUCL_CLASSIFY_SE_LOG,
    threads: KRAKEN2_CLASSIFY_THREADS
    wrapper:
        KRAKEN2_CLASSIFY_WRAPPER


rule kraken2_prot_read_classif_pe:
    input:
        fqs=[
            KRAKEN2_NUCL_UNCLASSIF_FASTQ_R1_FILE,
            KRAKEN2_NUCL_UNCLASSIF_FASTQ_R2_FILE,
        ],
        db=KRAKEN2_PROT_DB_DIR,
        build_done=KRAKEN2_PROT_DB_BUILD_DONE_FILE,
    params:
        output="-",
        paired_end=True,
        extra=(
            f"{config['kraken2']['classify']['extra']['paired_end']} "
            f"{config['kraken2']['classify']['extra']['common']}"
        ),
    output:
        classif=[
            temp(KRAKEN2_PROT_CLASSIF_FASTQ_R1_DC_FILE),
            temp(KRAKEN2_PROT_CLASSIF_FASTQ_R2_DC_FILE),
        ],
        unclassif=[
            temp(KRAKEN2_PROT_UNCLASSIF_FASTQ_R1_DC_FILE),
            temp(KRAKEN2_PROT_UNCLASSIF_FASTQ_R2_DC_FILE),
        ],
        report=KRAKEN2_PROT_REPORT_PE_FILE,
    log:
        KRAKEN2_PROT_CLASSIFY_PE_LOG,
    threads: KRAKEN2_CLASSIFY_THREADS
    wrapper:
        KRAKEN2_CLASSIFY_WRAPPER


rule kraken2_prot_read_classif_se:
    input:
        fqs=KRAKEN2_NUCL_UNCLASSIF_FASTQ_SE_FILE,
        db=KRAKEN2_PROT_DB_DIR,
        build_done=KRAKEN2_PROT_DB_BUILD_DONE_FILE,
    params:
        output="-",
        paired_end=False,
        extra=config["kraken2"]["classify"]["extra"]["common"],
    output:
        classif=temp(KRAKEN2_PROT_CLASSIF_FASTQ_SE_DC_FILE),
        unclassif=temp(KRAKEN2_PROT_UNCLASSIF_FASTQ_SE_DC_FILE),
        report=KRAKEN2_PROT_REPORT_SE_FILE,
    log:
        KRAKEN2_PROT_CLASSIFY_SE_LOG,
    threads: KRAKEN2_CLASSIFY_THREADS
    wrapper:
        KRAKEN2_CLASSIFY_WRAPPER


rule kraken2_compressed_fastq:
    input:
        KRAKEN2_FASTQ_DC_FILE,
    output:
        KRAKEN2_FASTQ_CP_FILE,
    log:
        KRAKEN2_COMPRESSED_FASTQ_LOG,
    threads: PIGZ_THREADS
    wrapper:
        PIGZ_WRAPPER


rule kraken2_combined_report:
    input:
        lambda wc: (
            [KRAKEN2_NUCL_REPORT_PE_FILE, KRAKEN2_PROT_REPORT_PE_FILE]
            if wc.etype == "pe"
            else [KRAKEN2_NUCL_REPORT_SE_FILE, KRAKEN2_PROT_REPORT_SE_FILE]
        ),
    params:
        extra=config["krakentools"]["combine_kreports"]["extra"],
    output:
        KRAKEN2_COMBINED_REPORT_FILE,
    log:
        KRAKEN2_COMBINED_REPORT_LOG,
    wrapper:
        KRAKENTOOLS_COMBINE_KREPORTS_WRAPPER


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
        task="build",
        extra=config["krakenuniq"]["build"]["extra"],
    output:
        touch(KRAKENUNIQ_DB_BUILD_DONE_FILE),
    log:
        KRAKENUNIQ_DB_BUILD_LOG,
    threads: KRAKENUNIQ_BUILD_THREADS
    wrapper:
        KRAKENUNIQ_BUILD_WRAPPER


rule krakenuniq_read_classif_pe:
    input:
        fqs=[BOWTIE2_FILTERED_FASTQ_R1_FILE, BOWTIE2_FILTERED_FASTQ_R2_FILE],
        db=KRAKENUNIQ_DB_DIR,
        build_done=KRAKENUNIQ_DB_BUILD_DONE_FILE,
    params:
        output="off",
        paired_end=True,
        extra=(
            f"{config['krakenuniq']['classify']['extra']['paired_end']} "
            f"{config['krakenuniq']['classify']['extra']['common']}"
        ),
    output:
        classif=[
            KRAKENUNIQ_CLASSIF_FASTQ_R1_FILE,
            KRAKENUNIQ_CLASSIF_FASTQ_R2_FILE,
        ],
        unclassif=[
            KRAKENUNIQ_UNCLASSIF_FASTQ_R1_FILE,
            KRAKENUNIQ_UNCLASSIF_FASTQ_R2_FILE,
        ],
        report=KRAKENUNIQ_REPORT_PE_FILE,
    log:
        KRAKENUNIQ_CLASSIFY_PE_LOG,
    threads: KRAKENUNIQ_CLASSIFY_THREADS
    wrapper:
        KRAKENUNIQ_CLASSIFY_WRAPPER


rule krakenuniq_read_classif_se:
    input:
        fqs=BOWTIE2_FILTERED_FASTQ_SE_FILE,
        db=KRAKENUNIQ_DB_DIR,
        build_done=KRAKENUNIQ_DB_BUILD_DONE_FILE,
    params:
        output="off",
        paired_end=False,
        extra=config["krakenuniq"]["classify"]["extra"]["common"],
    output:
        classif=KRAKENUNIQ_CLASSIF_FASTQ_SE_FILE,
        unclassif=KRAKENUNIQ_UNCLASSIF_FASTQ_SE_FILE,
        report=KRAKENUNIQ_REPORT_SE_FILE,
    log:
        KRAKENUNIQ_CLASSIFY_SE_LOG,
    threads: KRAKENUNIQ_CLASSIFY_THREADS
    wrapper:
        KRAKENUNIQ_CLASSIFY_WRAPPER
