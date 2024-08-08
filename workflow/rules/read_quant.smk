rule bracken_read_quant:
    input:
        report=lambda wc: (
            KRAKENUNIQ_REPORT_FILE
            if KRAKEN_MODE == "krakenuniq"
            else (
                KRAKEN2_COMBINED_REPORT_FILE
                if KRAKEN_MODE == "kraken2" and KRAKEN2_TSEARCH_UNCLASSIF
                else KRAKEN2_NUCL_REPORT_FILE
            )
        ),
        kdb_done=(
            KRAKEN2_NUCL_DB_DONE_FILE
            if KRAKEN_MODE == "kraken2"
            else KRAKENUNIQ_DB_DONE_FILE
        ),
        bdb_done=expand(BRACKEN_DB_DONE_FILE, **EXPAND_PARAMS),
        # readlen=READ_LENGTH_FILE,
    params:
        db=KRAKEN2_NUCL_DB_DIR if KRAKEN_MODE == "kraken2" else KRAKENUNIQ_DB_DIR,
        readlen=lambda wc: int(
            GDC_READGRP_META_DF.loc[wc.rg_id, "read_length"]
            if GDC_BAM_META_DF.loc[wc.bam_id, "num_uniq_read_groups"] > 1
            else GDC_BAM_META_DF.loc[wc.bam_id, "read_length"]
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


def bracken_count_files(wc):
    rg_ids, sfxs = glob_wildcards(
        join(
            GDC_FASTQ_RESULTS_DIR,
            wc.bam_id,
            "{rg_id,[0-9a-f\-]{36}}_unmapped_{sfx,(1|2|s){1}}.fq.gz",
        )
    )
    return expand(
        join(BRACKEN_QUANT_RESULTS_DIR, wc.bam_id, "{rg_id}_counts.tsv"),
        rg_id=rg_ids,
    )


rule bracken_combined_counts:
    input:
        bracken_count_files,
    output:
        BRACKEN_COMBINED_COUNT_FILE,
    log:
        BRACKEN_COMBINED_COUNT_LOG,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/braken_combined_counts.py"


rule bracken_count_matrix:
    input:
        expand(BRACKEN_COMBINED_COUNT_FILE, **EXPAND_PARAMS),
    params:
        samples=EXPAND_PARAMS["bam_id"],
    output:
        BRACKEN_COUNT_MATRIX_FILE,
    log:
        BRACKEN_COUNT_MATRIX_LOG,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/bracken_count_matrix.py"
