rule samtools_fastq:
    input:
        INPUT_BAM_FILE,
    params:
        extra=config["samtools"]["fastq"]["extra"],
    output:
        INPUT_FASTQ_FILE,
    log:
        INPUT_FASTQ_LOG,
    threads: SAMTOOLS_FASTQ_THREADS
    wrapper:
        SAMTOOLS_FASTQ_WRAPPER
