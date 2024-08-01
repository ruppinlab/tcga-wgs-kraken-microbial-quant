rule gdc_unmapped_bam:
    params:
        url="https://api.gdc.cancer.gov/slicing/view/{bam_id}?region=unmapped",
        headers={"X-Auth-Token": GDC_TOKEN},
        method="GET",
    output:
        temp(GDC_UNMAPPED_BAM_FILE),
    log:
        GDC_UNMAPPED_BAM_LOG,
    message:
        "{params.url}"
    retries: config["download"]["retries"]
    script:
        "../scripts/url_file.py"


rule gdc_unmapped_fastq_pe:
    input:
        GDC_UNMAPPED_BAM_FILE,
    params:
        outdir=lambda wc: GDC_FILE_RESULTS_DIR if wc.level == "rg" else None,
        paired_end=True,
        per_readgrp=lambda wc: True if wc.level == "rg" else False,
        O=GDC_UNMAPPED_FASTQ_O1_FILE,
        O2=GDC_UNMAPPED_FASTQ_O2_FILE,
        extra=lambda wc: (
            (
                f"{config['biobambam2']['bamtofastq']['extra']['common']} "
                f"{config['biobambam2']['bamtofastq']['extra']['paired_end']} "
                f"{config['biobambam2']['bamtofastq']['extra']['per_readgrp']}"
            )
            if wc.level == "rg"
            else (
                f"{config['biobambam2']['bamtofastq']['extra']['common']} "
                f"{config['biobambam2']['bamtofastq']['extra']['paired_end']}"
            )
        ),
    output:
        F=temp(GDC_UNMAPPED_FASTQ_R1_FILE),
        F2=temp(GDC_UNMAPPED_FASTQ_R2_FILE),
    log:
        GDC_UNMAPPED_FASTQ_LOG,
    wrapper:
        BIOBAMBAM2_BAMTOFASTQ_WRAPPER


rule gdc_unmapped_fastq_se:
    input:
        GDC_UNMAPPED_BAM_FILE,
    params:
        outdir=lambda wc: GDC_FILE_RESULTS_DIR if wc.level == "rg" else None,
        paired_end=False,
        per_readgrp=lambda wc: True if wc.level == "rg" else False,
        extra=lambda wc: (
            (
                f"{config['biobambam2']['bamtofastq']['extra']['common']} "
                f"{config['biobambam2']['bamtofastq']['extra']['per_readgrp']}"
            )
            if wc.level == "rg"
            else config["biobambam2"]["bamtofastq"]["extra"]["common"]
        ),
    output:
        S=temp(GDC_UNMAPPED_FASTQ_SE_FILE),
    log:
        GDC_UNMAPPED_FASTQ_LOG,
    wrapper:
        BIOBAMBAM2_BAMTOFASTQ_WRAPPER
