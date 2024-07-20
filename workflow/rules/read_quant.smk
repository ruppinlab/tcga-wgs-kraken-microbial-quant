rule bracken_db_build:
    input:
        KRAKEN2_DB_DIR if KRAKEN_MODE == "kraken2" else KRAKENUNIQ_DB_DIR,
    params:
        klen=35 if KRAKEN_MODE == "kraken2" else 31,
        ktype=KRAKEN_MODE,
        readlen="{readlen}",
    output:
        touch(BRACKEN_DB_BUILD_DONE_FILE),
    log:
        BRACKEN_DB_BUILD_LOG,
    threads: BRACKEN_BUILD_THREADS
    wrapper:
        BRACKEN_BUILD_WRAPPER


rule bracken_read_quant:
    input:
        db=KRAKEN2_DB_DIR if KRAKEN_MODE == "kraken2" else KRAKENUNIQ_DB_DIR,
        report=(
            KRAKEN2_REPORT_FILE if KRAKEN_MODE == "kraken2" else KRAKENUNIQ_REPORT_FILE
        ),
        build_done=expand(BRACKEN_DB_BUILD_DONE_FILE, **EXPAND_PARAMS),
    params:
        readlen=lambda wc: int(
            GDC_BAM_META_DF.loc[wc.bam_id, "read_length"]
            if wc.method == "sg"
            else GDC_READGRP_META_DF.loc[wc.rg_id, "read_length"]
        ),
        db_readlens=BRACKEN_DB_READ_LENGTHS,
        level=config["bracken"]["quant"]["level"],
        threshold=config["bracken"]["quant"]["threshold"],
    output:
        counts=BRACKEN_COUNT_FILE,
        report=BRACKEN_REPORT_FILE,
    log:
        BRACKEN_QUANT_LOG,
    wrapper:
        BRACKEN_QUANT_WRAPPER


rule bracken_merged_rg_counts:
    conda:
        "../envs/pandas.yaml"
    input:
        lambda wc: GDC_READGRP_META_DF[GDC_READGRP_META_DF["file_id"] == wc.bam_id]
        .apply(
            lambda x: join(
                BRACKEN_QUANT_RESULTS_DIR,
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
        BRACKEN_BAM_COUNT_FILES,
    params:
        samples=BRACKEN_BAM_ALIQUOT_IDS,
    output:
        BRACKEN_COUNT_MATRIX_FILE,
    log:
        BRACKEN_COUNT_MATRIX_LOG,
    script:
        "../scripts/braken_count_matrix.py"
