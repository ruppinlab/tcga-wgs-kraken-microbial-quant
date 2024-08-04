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
    conda:
        "../envs/samtools.yaml"
    script:
        "../scripts/url_file.py"


rule gdc_unmapped_fastq_pe:
    input:
        GDC_UNMAPPED_BAM_FILE,
    params:
        outputdir=lambda wc: (
            join(GDC_RESULTS_DIR, wc.bam_id, wc.level) if wc.level == "rg" else None
        ),
        paired_end=True,
        per_readgrp=lambda wc: True if wc.level == "rg" else False,
        rg_meta_df=lambda wc: (
            GDC_READGRP_META_DF.loc[GDC_READGRP_META_DF["file_id"] == wc.bam_id]
            if wc.level == "rg"
            else None
        ),
        O=GDC_UNMAPPED_FASTQ_O1_FILE,
        O2=GDC_UNMAPPED_FASTQ_O2_FILE,
        outputperreadgroupsuffixF="_unmapped_1.fq.gz",
        outputperreadgroupsuffixF2="_unmapped_2.fq.gz",
        outputperreadgroupsuffixO="_unmapped_o1.fq.gz",
        outputperreadgroupsuffixO2="_unmapped_o2.fq.gz",
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
        outputdir=lambda wc: (
            join(GDC_RESULTS_DIR, wc.bam_id, wc.level) if wc.level == "rg" else None
        ),
        paired_end=False,
        per_readgrp=lambda wc: True if wc.level == "rg" else False,
        rg_meta_df=lambda wc: (
            GDC_READGRP_META_DF.loc[GDC_READGRP_META_DF["file_id"] == wc.bam_id]
            if wc.level == "rg"
            else None
        ),
        outputperreadgroupsuffixS="_unmapped_s.fq.gz",
        extra=lambda wc: (
            (
                f"{config['biobambam2']['bamtofastq']['extra']['common']} "
                f"{config['biobambam2']['bamtofastq']['extra']['per_readgrp']}"
            )
            if wc.level == "rg"
            else config["biobambam2"]["bamtofastq"]["extra"]["common"]
        ),
    output:
        temp(GDC_UNMAPPED_FASTQ_SE_FILE),
    log:
        GDC_UNMAPPED_FASTQ_LOG,
    wrapper:
        BIOBAMBAM2_BAMTOFASTQ_WRAPPER
