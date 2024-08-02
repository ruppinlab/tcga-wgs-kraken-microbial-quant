rule bbmap_read_length_histogram:
    input:
        lambda wc: (
            [HOST_FILTERED_FASTQ_R1_FILE, HOST_FILTERED_FASTQ_R2_FILE]
            if wc.etype == "pe"
            else HOST_FILTERED_FASTQ_SE_FILE
        ),
    params:
        extra=config["bbmap"]["readlength"]["extra"],
    output:
        READ_LENGTH_HISTOGRAM_FILE,
    log:
        READ_LENGTH_HISTOGRAM_LOG,
    threads: BBMAP_READLENGTH_THREADS
    wrapper:
        BBMAP_READLENGTH_WRAPPER


rule bbmap_max_read_length:
    conda:
        "../envs/pandas.yaml"
    input:
        READ_LENGTH_HISTOGRAM_FILE,
    output:
        READ_LENGTH_FILE,
    log:
        READ_LENGTH_LOG,
    script:
        "../scripts/read_length.py"
