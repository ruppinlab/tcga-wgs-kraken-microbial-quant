rule bowtie2_chm13_index:
    input:
        ref=CHM13_FASTA_FILE,
    output:
        idx=multiext(
            BOWTIE2_INDEX_DIR,
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        BOWTIE2_INDEX_LOG,
    threads: BOWTIE2_BUILD_THREADS
    wrapper:
        BOWTIE2_BUILD_WRAPPER


rule bowtie2_chm13_align:
    input:
        sample=INPUT_FASTQ_FILE,
        idx=multiext(
            BOWTIE2_INDEX_DIR,
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        unaligned=BOWTIE2_UNMAPPED_SAM_FILE,
    log:
        BOWTIE2_SAM_LOG,
    threads: BOWTIE2_ALIGN_THREADS
    wrapper:
        BOWTIE2_ALIGN_WRAPPER
