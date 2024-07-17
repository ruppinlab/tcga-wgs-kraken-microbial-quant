rule krakenuniq_classify_microbialdb_pe:
    input:
        fqs=[BOWTIE2_FILTERED_FASTQ_R1_FILE, BOWTIE2_FILTERED_FASTQ_R2_FILE],
        db=config["resources"]["krakenuniq"]["microbialdb"]["dir"],
    params:
        extra=lambda wc: (
            f"--paired --check-names config['krakenuniq']['extra']"
            if GDC_READGRP_META_DF.loc[wc.rg_id, "is_paired_end"]
            else config["krakenuniq"]["extra"]
        ),
    output:
        KRAKENUNIQ_MICROBIALDB_CLASSIF_PE_FILE,
        report=KRAKENUNIQ_MICROBIALDB_REPORT_PE_FILE,
    log:
        KRAKENUNIQ_MICROBIALDB_CLASSIF_LOG,
    threads: KRAKENUNIQ_THREADS
    wrapper:
        KRAKENUNIQ_WRAPPER


rule krakenuniq_classify_microbialdb_se:
    input:
        fqs=BOWTIE2_FILTERED_FASTQ_SE_FILE,
        db=config["resources"]["krakenuniq"]["microbialdb"]["dir"],
    params:
        extra=lambda wc: (
            f"--paired --check-names config['krakenuniq']['extra']"
            if GDC_READGRP_META_DF.loc[wc.rg_id, "is_paired_end"]
            else config["krakenuniq"]["extra"]
        ),
    output:
        KRAKENUNIQ_MICROBIALDB_CLASSIF_SE_FILE,
        report=KRAKENUNIQ_MICROBIALDB_REPORT_SE_FILE,
    log:
        KRAKENUNIQ_MICROBIALDB_CLASSIF_LOG,
    threads: KRAKENUNIQ_THREADS
    wrapper:
        KRAKENUNIQ_WRAPPER


rule bracken_quantify_microbialdb_pe:
    input:
        db=config["resources"]["krakenuniq"]["microbialdb"]["dir"],
        report=KRAKENUNIQ_MICROBIALDB_REPORT_PE_FILE,
        readlen=READ_LENGTH_FILE,
    params:
        # readlen=lambda wc: GDC_READGRP_META_DF.loc[wc.rg_id, "read_length"],
        extra=config["bracken"]["extra"],
    output:
        BRACKEN_MICROBIALDB_QUANT_PE_FILE,
    log:
        BRACKEN_MICROBIALDB_QUANT_LOG,
    wrapper:
        BRACKEN_WRAPPER


rule bracken_quantify_microbialdb_se:
    input:
        db=config["resources"]["krakenuniq"]["microbialdb"]["dir"],
        report=KRAKENUNIQ_MICROBIALDB_REPORT_SE_FILE,
        readlen=READ_LENGTH_FILE,
    params:
        # readlen=lambda wc: GDC_READGRP_META_DF.loc[wc.rg_id, "read_length"],
        extra=config["bracken"]["extra"],
    output:
        BRACKEN_MICROBIALDB_QUANT_SE_FILE,
    log:
        BRACKEN_MICROBIALDB_QUANT_LOG,
    wrapper:
        BRACKEN_WRAPPER
