rule krakenuniq_classify_microbialdb:
    input:
        fq1=BOWTIE2_FILTERED_FASTQ1_FILE,
        fq2=BOWTIE2_FILTERED_FASTQ2_FILE,
        db=config["resources"]["krakenuniq"]["microbialdb"]["dir"],
    params:
        extra=config["krakenuniq"]["extra"],
    output:
        KRAKENUNIQ_MICROBIALDB_CLASSIF_FILE,
        report=KRAKENUNIQ_MICROBIALDB_REPORT_FILE,
    log:
        KRAKENUNIQ_MICROBIALDB_CLASSIF_LOG,
    threads: KRAKENUNIQ_THREADS
    wrapper:
        KRAKENUNIQ_WRAPPER


rule bracken_quantify_microbialdb:
    input:
        db=config["resources"]["krakenuniq"]["microbialdb"]["dir"],
        report=KRAKENUNIQ_MICROBIALDB_REPORT_FILE,
    params:
        # readlen=,
        level=config["bracken"]["level"],
        threshold=config["bracken"]["threshold"],
        extra=config["bracken"]["extra"],
    output:
        BRACKEN_MICROBIALDB_QUANT_FILE,
    log:
        BRACKEN_MICROBIALDB_QUANT_LOG,
    wrapper:
        BRACKEN_WRAPPER
