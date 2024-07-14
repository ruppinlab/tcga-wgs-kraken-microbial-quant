rule gdc_unmapped_bam_file:
    conda:
        "../envs/httr.yaml"
    params:
        token=config["input"]["gdc"]["token"],
        bam_id=GDC_BAM_ID_WILDCARD_STR,
    output:
        GDC_UNMAPPED_BAM_FILE,
    log:
        GDC_UNMAPPED_BAM_LOG,
    retries: config["download"]["retries"]
    script:
        "../scripts/gdc_unmapped_bam.R"


rule gdc_unmapped_fastq_pe:
    input:
        lambda wc: f"{GDC_READGRP_META_DF.loc[wc.rg_id].file_id}_unmapped.bam",
    params:
        outdir=GDC_RESULTS_DIR,
        extra=config["biobambam2"]["bamtofastq"]["extra"],
    output:
        GDC_UNMAPPED_FASTQ_1_FILE,
        GDC_UNMAPPED_FASTQ_2_FILE,
    log:
        GDC_UNMAPPED_FASTQ_LOG,
    wrapper:
        BIOBAMBAM_BAMTOFASTQ_WRAPPER


rule gdc_unmapped_fastq_se:
    input:
        lambda wc: f"{GDC_READGRP_META_DF.loc[wc.rg_id].file_id}_unmapped.bam",
    params:
        outdir=GDC_RESULTS_DIR,
        extra=config["biobambam2"]["bamtofastq"]["extra"],
    output:
        GDC_UNMAPPED_FASTQ_S_FILE,
    log:
        GDC_UNMAPPED_FASTQ_LOG,
    wrapper:
        BIOBAMBAM_BAMTOFASTQ_WRAPPER
