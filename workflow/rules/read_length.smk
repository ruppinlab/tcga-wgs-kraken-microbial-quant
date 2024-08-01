rule bbmap_read_length_histogram:
    input:
        lambda wc: (
            [HOST_FILTERED_FASTQ_R1_FILE, HOST_FILTERED_FASTQ_R2_FILE]
            if wc.etype == "pe"
            else HOST_FILTERED_FASTQ_SE_FILE
        ),
    output:
        READ_LENGTH_HISTOGRAM_FILE,
    log:
        READ_LENGTH_HISTOGRAM_LOG,
    wrapper:
        READ_LENGTH_WRAPPER


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
