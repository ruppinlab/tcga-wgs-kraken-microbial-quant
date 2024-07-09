rule gdc_bam_meta:
    conda:
        "../envs/httr.yaml"
    params:
        program_names=config["input"]["gdc"]["program_names"],
        sample_types=config["input"]["gdc"]["sample_types"],
        exp_strategy=config["input"]["gdc"]["exp_strategy"],
        workflow_types=config["input"]["gdc"]["workflow_types"],
    output:
        GDC_BAM_META_FILE,
    log:
        GDC_BAM_META_LOG,
    script:
        "../scripts/gdc_file_meta.py"


checkpoint gdc_file_touch:
    conda:
        "../envs/pandas.yaml"
    input:
        GDC_BAM_META_FILE,
    output:
        temp(directory(GDC_BAM_TOUCH_DIR)),
    log:
        GDC_BAM_TOUCH_LOG,
    script:
        "../scripts/gdc_file_touch.py"


rule gdc_unmapped_bam_file:
    conda:
        "../envs/httr.yaml"
    input:
        GDC_BAM_TOUCH_FILE,
    output:
        GDC_BAM_FILE,
    log:
        GDC_BAM_LOG,
    retries: config["download"]["retries"]
    script:
        "../scripts/gdc_unmapped_bam.R"


rule gdc_sorted_bam:
    input:
        GDC_BAM_FILE,
    params:
        extra=config["samtools"]["sort"]["extra"],
    output:
        GDC_SORTED_BAM_FILE,
    log:
        GDC_SORTED_BAM_LOG,
    threads: SAMTOOLS_SORT_THREADS
    wrapper:
        SAMTOOLS_SORT_WRAPPER


rule gdc_fastq:
    input:
        GDC_SORTED_BAM_FILE,
    params:
        extra=config["samtools"]["fastq"]["extra"],
    output:
        GDC_FASTQ_FILE,
    log:
        GDC_FASTQ_LOG,
    threads: SAMTOOLS_FASTQ_THREADS
    wrapper:
        SAMTOOLS_FASTQ_WRAPPER
