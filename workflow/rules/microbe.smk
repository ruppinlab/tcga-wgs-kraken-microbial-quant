rule krakenuniq_classify:
    input:
        fqs=lambda wc: (
            [BOWTIE2_FILTERED_FASTQ_R1_FILE, BOWTIE2_FILTERED_FASTQ_R2_FILE]
            if wc.etype == "pe"
            else BOWTIE2_FILTERED_FASTQ_SE_FILE
        ),
        db=config["resources"]["krakenuniq"]["microbialdb"]["dir"],
    params:
        extra=lambda wc: (
            f"--paired --check-names {config['krakenuniq']['extra']}"
            if wc.etype == "pe"
            else config["krakenuniq"]["extra"]
        ),
    output:
        KRAKENUNIQ_CLASSIF_FILE,
        report=KRAKENUNIQ_REPORT_FILE,
    log:
        KRAKENUNIQ_CLASSIF_LOG,
    threads: KRAKENUNIQ_THREADS
    wrapper:
        KRAKENUNIQ_WRAPPER


rule bracken_quantify:
    input:
        db=config["resources"]["krakenuniq"]["microbialdb"]["dir"],
        report=KRAKENUNIQ_REPORT_FILE,
        readlen=READ_LENGTH_FILE,
    params:
        # readlen=lambda wc: GDC_READGRP_META_DF.loc[wc.rg_id, "read_length"],
        extra=config["bracken"]["extra"],
    output:
        BRACKEN_QUANT_FILE,
    log:
        BRACKEN_QUANT_LOG,
    wrapper:
        BRACKEN_WRAPPER
