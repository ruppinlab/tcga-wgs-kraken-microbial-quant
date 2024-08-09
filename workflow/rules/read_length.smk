rule bbmap_read_length_histogram_pe:
    input:
        [HOST_FILTERED_FASTQ_R1_FILE, HOST_FILTERED_FASTQ_R2_FILE],
    params:
        extra=config["bbmap"]["readlength"]["extra"],
    output:
        READ_LENGTH_HISTOGRAM_FILE,
    log:
        READ_LENGTH_HISTOGRAM_LOG,
    group:
        "group_{bam_id}"
    threads: BBMAP_READLENGTH_THREADS
    wrapper:
        BBMAP_READLENGTH_WRAPPER


rule bbmap_read_length_histogram_se:
    input:
        HOST_FILTERED_FASTQ_SE_FILE,
    params:
        extra=config["bbmap"]["readlength"]["extra"],
    output:
        READ_LENGTH_HISTOGRAM_FILE,
    log:
        READ_LENGTH_HISTOGRAM_LOG,
    group:
        "group_{bam_id}"
    threads: BBMAP_READLENGTH_THREADS
    wrapper:
        BBMAP_READLENGTH_WRAPPER


rule bbmap_max_read_length:
    input:
        READ_LENGTH_HISTOGRAM_FILE,
    output:
        READ_LENGTH_FILE,
    log:
        READ_LENGTH_LOG,
    group:
        "group_{bam_id}"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/read_length.py"
