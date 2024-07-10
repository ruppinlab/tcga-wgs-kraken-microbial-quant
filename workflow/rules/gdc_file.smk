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


rule gdc_sorted_unmapped_bam:
    input:
        GDC_UNMAPPED_BAM_FILE,
    params:
        extra=config["samtools"]["sort"]["extra"],
    output:
        GDC_SORTED_UNMAPPED_BAM_FILE,
    log:
        GDC_SORTED_UNMAPPED_BAM_LOG,
    threads: SAMTOOLS_SORT_THREADS
    wrapper:
        SAMTOOLS_SORT_WRAPPER


rule gdc_unmapped_fastq:
    input:
        GDC_SORTED_UNMAPPED_BAM_FILE,
    params:
        extra=config["samtools"]["fastq"]["extra"],
    output:
        GDC_UNMAPPED_FASTQ1_FILE,
        GDC_UNMAPPED_FASTQ2_FILE,
    log:
        GDC_UNMAPPED_FASTQ_LOG,
    threads: SAMTOOLS_FASTQ_THREADS
    wrapper:
        SAMTOOLS_FASTQ_SEPARATE_WRAPPER
