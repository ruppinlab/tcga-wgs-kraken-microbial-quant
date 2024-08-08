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
    resources:
        gdc_download_jobs=1,
    retries: config["download"]["retries"]
    conda:
        "../envs/samtools.yaml"
    script:
        "../scripts/url_bam_file.py"


checkpoint gdc_unmapped_fastqs:
    input:
        GDC_UNMAPPED_BAM_FILE,
    params:
        rg_meta_df=lambda wc: GDC_READGRP_META_DF.loc[
            GDC_READGRP_META_DF["file_id"] == wc.bam_id
        ],
        suffixes={
            "F": "_unmapped_1.fq.gz",
            "F2": "_unmapped_2.fq.gz",
            "O": "_unmapped_o1.fq.gz",
            "O2": "_unmapped_o2.fq.gz",
            "S": "_unmapped_s.fq.gz",
        },
        extra=config["biobambam2"]["bamtofastq"]["extra"],
    output:
        directory(GDC_UNMAPPED_FASTQ_DIR),
    log:
        GDC_UNMAPPED_FASTQ_LOG,
    wrapper:
        BIOBAMBAM2_BAMTOFASTQ_WRAPPER
