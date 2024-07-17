rule gdc_unmapped_bam:
    params:
        token=config["input"]["gdc"]["token"],
        bam_id=GDC_BAM_ID_WILDCARD_STR,
    output:
        GDC_UNMAPPED_BAM_FILE,
    log:
        GDC_UNMAPPED_BAM_LOG,
    retries: config["download"]["retries"]
    script:
        "../scripts/gdc_unmapped_bam.py"


rule gdc_unmapped_fastq_pe:
    input:
        GDC_UNMAPPED_BAM_FILE,
    params:
        paired_end=True,
        extra=(
            f"{config['biobambam2']['bamtofastq']['extra']['common']} "
            f"{config['biobambam2']['bamtofastq']['extra']['paired_end']}"
        ),
    output:
        GDC_UNMAPPED_FASTQ_R1_FILE,
        GDC_UNMAPPED_FASTQ_R2_FILE,
        temp(GDC_UNMAPPED_FASTQ_O1_FILE),
        temp(GDC_UNMAPPED_FASTQ_O2_FILE),
    log:
        GDC_UNMAPPED_FASTQ_LOG,
    wrapper:
        BIOBAMBAM2_BAMTOFASTQ_WRAPPER


rule gdc_unmapped_fastq_pe_per_readgrp:
    input:
        GDC_UNMAPPED_BAM_FILE,
    params:
        outdir=GDC_BAM_RESULTS_DIR,
        paired_end=True,
        per_readgrp=True,
        extra=(
            f"{config['biobambam2']['bamtofastq']['extra']['common']} "
            f"{config['biobambam2']['bamtofastq']['extra']['paired_end']}"
            f"{config['biobambam2']['bamtofastq']['extra']['per_readgrp']}"
        ),
    output:
        GDC_UNMAPPED_FASTQ_R1_FILE,
        GDC_UNMAPPED_FASTQ_R2_FILE,
        temp(GDC_UNMAPPED_FASTQ_O1_FILE),
        temp(GDC_UNMAPPED_FASTQ_O2_FILE),
    log:
        GDC_UNMAPPED_FASTQ_LOG,
    wrapper:
        BIOBAMBAM2_BAMTOFASTQ_WRAPPER


rule gdc_unmapped_fastq_se:
    input:
        GDC_UNMAPPED_BAM_FILE,
    params:
        extra=config["biobambam2"]["bamtofastq"]["extra"]["common"],
    output:
        GDC_UNMAPPED_FASTQ_SE_FILE,
    log:
        GDC_UNMAPPED_FASTQ_LOG,
    wrapper:
        BIOBAMBAM2_BAMTOFASTQ_WRAPPER


rule gdc_unmapped_fastq_se_per_readgrp:
    input:
        GDC_UNMAPPED_BAM_FILE,
    params:
        outdir=GDC_BAM_RESULTS_DIR,
        per_readgrp=True,
        extra=(
            f"{config['biobambam2']['bamtofastq']['extra']['common']} "
            f"{config['biobambam2']['bamtofastq']['extra']['per_readgrp']}"
        ),
    output:
        GDC_UNMAPPED_FASTQ_SE_FILE,
    log:
        GDC_UNMAPPED_FASTQ_LOG,
    wrapper:
        BIOBAMBAM2_BAMTOFASTQ_WRAPPER
