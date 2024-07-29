rule kraken2_nucl_read_classif_pe:
    input:
        fqs=[BOWTIE2_FILTERED_FASTQ_R1_FILE, BOWTIE2_FILTERED_FASTQ_R2_FILE],
        db_done=KRAKEN2_NUCL_DB_DONE_FILE,
    params:
        db=KRAKEN2_NUCL_DB_DIR,
        output="-",
        paired_end=True,
        extra=(
            f"{config['kraken2']['classify']['extra']['paired_end']} "
            f"{config['kraken2']['classify']['extra']['common']}"
        ),
    output:
        classif=[
            temp(KRAKEN2_NUCL_CLASSIF_FASTQ_R1_FILE),
            temp(KRAKEN2_NUCL_CLASSIF_FASTQ_R2_FILE),
        ],
        unclassif=[
            temp(KRAKEN2_NUCL_UNCLASSIF_FASTQ_R1_FILE),
            temp(KRAKEN2_NUCL_UNCLASSIF_FASTQ_R2_FILE),
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
        db_done=KRAKEN2_NUCL_DB_DONE_FILE,
    params:
        db=KRAKEN2_NUCL_DB_DIR,
        output="-",
        paired_end=False,
        extra=config["kraken2"]["classify"]["extra"]["common"],
    output:
        classif=temp(KRAKEN2_NUCL_CLASSIF_FASTQ_SE_FILE),
        unclassif=temp(KRAKEN2_NUCL_UNCLASSIF_FASTQ_SE_FILE),
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
        db_done=KRAKEN2_PROT_DB_DONE_FILE,
    params:
        db=KRAKEN2_PROT_DB_DIR,
        output="-",
        paired_end=True,
        extra=(
            f"{config['kraken2']['classify']['extra']['paired_end']} "
            f"{config['kraken2']['classify']['extra']['common']}"
        ),
    output:
        classif=[
            temp(KRAKEN2_PROT_CLASSIF_FASTQ_R1_FILE),
            temp(KRAKEN2_PROT_CLASSIF_FASTQ_R2_FILE),
        ],
        unclassif=[
            temp(KRAKEN2_PROT_UNCLASSIF_FASTQ_R1_FILE),
            temp(KRAKEN2_PROT_UNCLASSIF_FASTQ_R2_FILE),
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
        db_done=KRAKEN2_PROT_DB_DONE_FILE,
    params:
        db=KRAKEN2_PROT_DB_DIR,
        output="-",
        paired_end=False,
        extra=config["kraken2"]["classify"]["extra"]["common"],
    output:
        classif=temp(KRAKEN2_PROT_CLASSIF_FASTQ_SE_FILE),
        unclassif=temp(KRAKEN2_PROT_UNCLASSIF_FASTQ_SE_FILE),
        report=KRAKEN2_PROT_REPORT_SE_FILE,
    log:
        KRAKEN2_PROT_CLASSIFY_SE_LOG,
    threads: KRAKEN2_CLASSIFY_THREADS
    wrapper:
        KRAKEN2_CLASSIFY_WRAPPER


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


rule krakenuniq_read_classif_pe:
    input:
        fqs=[BOWTIE2_FILTERED_FASTQ_R1_FILE, BOWTIE2_FILTERED_FASTQ_R2_FILE],
        db_done=KRAKENUNIQ_DB_DONE_FILE,
    params:
        db=KRAKENUNIQ_DB_DIR,
        output="off",
        paired_end=True,
        extra=(
            f"{config['krakenuniq']['classify']['extra']['paired_end']} "
            f"{config['krakenuniq']['classify']['extra']['common']}"
        ),
    output:
        classif=[
            temp(KRAKENUNIQ_CLASSIF_FASTQ_R1_FILE),
            temp(KRAKENUNIQ_CLASSIF_FASTQ_R2_FILE),
        ],
        unclassif=[
            temp(KRAKENUNIQ_UNCLASSIF_FASTQ_R1_FILE),
            temp(KRAKENUNIQ_UNCLASSIF_FASTQ_R2_FILE),
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
        db_done=KRAKENUNIQ_DB_DONE_FILE,
    params:
        db=KRAKENUNIQ_DB_DIR,
        output="off",
        paired_end=False,
        extra=config["krakenuniq"]["classify"]["extra"]["common"],
    output:
        classif=temp(KRAKENUNIQ_CLASSIF_FASTQ_SE_FILE),
        unclassif=temp(KRAKENUNIQ_UNCLASSIF_FASTQ_SE_FILE),
        report=KRAKENUNIQ_REPORT_SE_FILE,
    log:
        KRAKENUNIQ_CLASSIFY_SE_LOG,
    threads: KRAKENUNIQ_CLASSIFY_THREADS
    wrapper:
        KRAKENUNIQ_CLASSIFY_WRAPPER
