rule samtools_fastq:
    input:
        INPUT_BAM_FILE,
    params:
        extra=config["samtools"]["fastq"]["extra"],
    output:
        BAM_FASTQ_FILE,
    log:
        BAM_FASTQ_LOG,
    threads: SAMTOOLS_THREADS
    wrapper:
        SAMTOOLS_WRAPPER
