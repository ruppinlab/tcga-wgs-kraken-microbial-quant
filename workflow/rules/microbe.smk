rule krakenuniq_classify_reads:
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


rule bracken_quantify_reads:
    input:
        db=config["resources"]["krakenuniq"]["microbialdb"]["dir"],
        report=KRAKENUNIQ_REPORT_FILE,
        readlen=READ_LENGTH_FILE,
    params:
        # readlen=lambda wc: GDC_READGRP_META_DF.loc[wc.rg_id, "read_length"],
        level=config["bracken"]["level"],
        threshold=config["bracken"]["threshold"],
    output:
        counts=BRACKEN_COUNT_FILE,
        report=BRACKEN_REPORT_FILE,
    log:
        BRACKEN_COUNT_LOG,
    wrapper:
        BRACKEN_WRAPPER


rule bracken_merged_rg_counts:
    conda:
        "../envs/pandas.yaml"
    input:
        lambda wc: GDC_READGRP_META_DF[GDC_READGRP_META_DF["file_id"] == wc.bam_id]
        .apply(
            lambda x: join(
                BRACKEN_COUNT_RESULTS_DIR,
                "rg",
                x["file_id"],
                "pe" if x["is_paired_end"] else "se",
                f"{x['read_group_id']}_counts.tsv",
            ),
            axis=1,
        )
        .tolist(),
    output:
        BRACKEN_MERGED_RG_COUNT_FILE,
    log:
        BRACKEN_MERGED_RG_COUNT_LOG,
    script:
        "../scripts/braken_sum_counts.py"


rule bracken_count_matrix:
    conda:
        "../envs/pandas.yaml"
    input:
        counts=BRACKEN_BAM_COUNT_FILES,
    params:
        samples=BRACKEN_BAM_IDS,
    output:
        BRACKEN_COUNT_MATRIX_FILE,
    log:
        BRACKEN_COUNT_MATRIX_LOG,
    script:
        "../scripts/braken_count_matrix.py"
