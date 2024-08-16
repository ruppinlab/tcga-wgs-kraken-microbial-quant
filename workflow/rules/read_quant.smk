rule bracken_read_quant:
    input:
        report=lambda wc: (
            KRAKEN2_COMBINED_REPORT_FILE
            if KRAKEN_MODE == "kraken2" and KRAKEN2_TSEARCH_UNCLASSIF
            else (
                KRAKEN2_NUCL_REPORT_FILE
                if KRAKEN_MODE == "kraken2"
                else KRAKENUNIQ_REPORT_FILE
            )
        ),
        kdb_done=(
            KRAKEN2_NUCL_DB_DONE_FILE
            if KRAKEN_MODE == "kraken2"
            else KRAKENUNIQ_DB_DONE_FILE
        ),
        bdb_done=expand(BRACKEN_DB_KMER_DISTR_DONE_FILE, **EXPAND_PARAMS),
        # readlen=READ_LENGTH_FILE,
    params:
        bracken=BRACKEN_QUANT_SCRIPT_PATH,
        db=KRAKEN2_NUCL_DB_DIR if KRAKEN_MODE == "kraken2" else KRAKENUNIQ_DB_DIR,
        readlen=lambda wc: int(
            GDC_READGRP_META_DF.loc[wc.rg_id, "read_length"]
            if GDC_BAM_META_DF.loc[wc.bam_id, "num_uniq_read_groups"] > 1
            else GDC_BAM_META_DF.loc[wc.bam_id, "read_length"]
        ),
        db_readlens=BRACKEN_DB_READ_LENGTHS,
        level=config["bracken"]["quant"]["level"],
        threshold=lambda wc: config["bracken"]["quant"]["threshold"][
            (
                "rg"
                if GDC_BAM_META_DF.loc[wc.bam_id, "num_uniq_read_groups"] > 1
                else "sg"
            )
        ],
    output:
        counts=BRACKEN_COUNT_FILE,
        report=BRACKEN_REPORT_FILE,
    log:
        BRACKEN_QUANT_LOG,
    # group:
    #     "{rg_id}"
    wrapper:
        BRACKEN_QUANT_WRAPPER


def bracken_rg_count_files(wildcards):
    gdc_rg_unmapped_fastq_dir = checkpoints.gdc_rg_unmapped_fastqs.get(
        **wildcards
    ).output[0]
    rg_ids, sfxs = glob_wildcards(
        join(
            gdc_rg_unmapped_fastq_dir,
            "{rg_id,[0-9a-f\\-]{36}}_unmapped_{sfx,(1|2|s){1}}.fq.gz",
        )
    )
    rg_id2etype = {}
    for rg_id, sfx in zip(rg_ids, sfxs):
        if rg_id not in rg_id2etype:
            rg_id2etype[rg_id] = {}
        rg_id2etype[rg_id]["pe" if sfx in ("1", "2") else "se"] = True
    rg_ids, etypes = [], []
    for k, v in rg_id2etype.items():
        for e in v.keys():
            rg_ids.append(k)
            etypes.append(e)
    return expand(
        join(
            BRACKEN_QUANT_RESULTS_DIR, wildcards.rg_bam_id, "{rg_id}_counts_{etype}.tsv"
        ),
        zip,
        rg_id=rg_ids,
        etype=etypes,
    )


rule bracken_combined_rg_counts:
    input:
        bracken_rg_count_files,
    output:
        BRACKEN_COMBINED_RG_COUNT_FILE,
    log:
        BRACKEN_COMBINED_RG_COUNT_LOG,
    # group:
    #     "{bam_id}"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/bracken_combined_counts.py"


rule bracken_count_matrix:
    input:
        BRACKEN_BAM_COUNT_FILES,
    params:
        samples=BRACKEN_BAM_IDS,
    output:
        BRACKEN_COUNT_MATRIX_FILE,
    log:
        BRACKEN_COUNT_MATRIX_LOG,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/bracken_count_matrix.py"
