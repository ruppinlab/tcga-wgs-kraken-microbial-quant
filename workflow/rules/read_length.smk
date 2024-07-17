rule bbmap_read_length_histogram_pe:
    input:
        fq1=BOWTIE2_FILTERED_FASTQ_R1_FILE,
        fq2=BOWTIE2_FILTERED_FASTQ_R2_FILE,
    output:
        READ_LENGTH_HISTOGRAM_PE_FILE,
    log:
        READ_LENGTH_HISTOGRAM_LOG,
    wrapper:
        READ_LENGTH_WRAPPER


rule bbmap_read_length_histogram_se:
    input:
        fq1=BOWTIE2_FILTERED_FASTQ_SE_FILE,
    output:
        READ_LENGTH_HISTOGRAM_SE_FILE,
    log:
        READ_LENGTH_HISTOGRAM_LOG,
    wrapper:
        READ_LENGTH_WRAPPER


rule bbmap_max_read_length_pe:
    conda:
        "../envs/pandas.yaml"
    input:
        READ_LENGTH_HISTOGRAM_PE_FILE,
    output:
        READ_LENGTH_PE_FILE,
    log:
        READ_LENGTH_LOG,
    script:
        "../scripts/read_length.py"


rule bbmap_max_read_length_se:
    conda:
        "../envs/pandas.yaml"
    input:
        READ_LENGTH_HISTOGRAM_SE_FILE,
    output:
        READ_LENGTH_SE_FILE,
    log:
        READ_LENGTH_LOG,
    script:
        "../scripts/read_length.py"
