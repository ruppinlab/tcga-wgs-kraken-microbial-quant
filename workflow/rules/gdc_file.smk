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


rule gdc_sg_unmapped_fastq_pe:
    input:
        GDC_SG_UNMAPPED_BAM_FILE,
    params:
        rg_meta_df=lambda wc: GDC_READGRP_META_DF.loc[
            GDC_READGRP_META_DF["file_id"] == wc.sg_bam_id
        ],
        O=GDC_SG_UNMAPPED_FASTQ_O1_FILE,
        O2=GDC_SG_UNMAPPED_FASTQ_O2_FILE,
        S=GDC_SG_UNMAPPED_FASTQ_SE_FILE,
        extra=config["biobambam2"]["bamtofastq"]["extra"],
    output:
        F=temp(GDC_SG_UNMAPPED_FASTQ_R1_FILE),
        F2=temp(GDC_SG_UNMAPPED_FASTQ_R2_FILE),
    log:
        GDC_SG_UNMAPPED_FASTQ_LOG,
    wrapper:
        BIOBAMBAM2_BAMTOFASTQ_WRAPPER


rule gdc_sg_unmapped_fastq_se:
    input:
        GDC_SG_UNMAPPED_BAM_FILE,
    params:
        rg_meta_df=lambda wc: GDC_READGRP_META_DF.loc[
            GDC_READGRP_META_DF["file_id"] == wc.sg_bam_id
        ],
        extra=config["biobambam2"]["bamtofastq"]["extra"],
    output:
        temp(GDC_SG_UNMAPPED_FASTQ_SE_FILE),
    log:
        GDC_SG_UNMAPPED_FASTQ_LOG,
    wrapper:
        BIOBAMBAM2_BAMTOFASTQ_WRAPPER


checkpoint gdc_rg_unmapped_fastqs:
    input:
        GDC_RG_UNMAPPED_BAM_FILE,
    params:
        rg_meta_df=lambda wc: GDC_READGRP_META_DF.loc[
            GDC_READGRP_META_DF["file_id"] == wc.rg_bam_id
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
        directory(GDC_RG_UNMAPPED_FASTQ_FILE_DIR),
    log:
        GDC_RG_UNMAPPED_FASTQ_LOG,
    wrapper:
        BIOBAMBAM2_BAMTOFASTQ_WRAPPER
