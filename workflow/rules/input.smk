rule samtools_unmapped_fastq:
    input:
        UNMAPPED_BAM_FILE,
    params:
        extra=config["samtools"]["unmapped"]["extra"],
    output:
        UNMAPPED_FASTQ_FILE,
    log:
        UNMAPPED_FASTQ_LOG,
    threads: SAMTOOLS_UNMAPPED_THREADS
    wrapper:
        SAMTOOLS_FASTQ_WRAPPER
