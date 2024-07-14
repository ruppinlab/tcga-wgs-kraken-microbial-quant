rule krakenuniq_classify_microbialdb:
    input:
        db=config["resources"]["krakenuniq"]["microbialdb"]["dir"],
    params:
        extra=lambda wc: (
            f"--paired --check-names config['krakenuniq']['extra']"
            if GDC_READGRP_META_DF.loc[wc.rg_id].is_paired_end
            else config["krakenuniq"]["extra"]
        ),
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
        readlen=lambda wc: GDC_READGRP_META_DF.loc[wc.rg_id].read_length,
        extra=config["bracken"]["extra"],
    output:
        BRACKEN_MICROBIALDB_QUANT_FILE,
    log:
        BRACKEN_MICROBIALDB_QUANT_LOG,
    wrapper:
        BRACKEN_WRAPPER
