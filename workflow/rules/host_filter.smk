rule bowtie2_index:
    input:
        ref=REF_FASTA_FILE,
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


rule bowtie2_align:
    input:
        sample=[GDC_UNMAPPED_FASTQ1_FILE, GDC_UNMAPPED_FASTQ2_FILE],
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
        temp(BOWTIE2_SAM_FILE),
        unaligned=BOWTIE2_FILTERED_SAM_FILE,
    log:
        BOWTIE2_ALIGN_LOG,
    threads: BOWTIE2_ALIGN_THREADS
    wrapper:
        BOWTIE2_ALIGN_WRAPPER


rule bowtie2_sorted_sam:
    input:
        BOWTIE2_FILTERED_SAM_FILE,
    params:
        extra=config["samtools"]["sort"]["extra"],
    output:
        BOWTIE2_SORTED_FILTERED_SAM_FILE,
    log:
        BOWTIE2_SORTED_FILTERED_SAM_LOG,
    threads: SAMTOOLS_SORT_THREADS
    wrapper:
        SAMTOOLS_SORT_WRAPPER


rule bowtie2_filtered_fastq:
    input:
        BOWTIE2_SORTED_FILTERED_SAM_FILE,
    params:
        extra=config["samtools"]["fastq"]["extra"],
    output:
        BOWTIE2_FILTERED_FASTQ1_FILE,
        BOWTIE2_FILTERED_FASTQ2_FILE,
    log:
        BOWTIE2_FILTERED_FASTQ_LOG,
    threads: SAMTOOLS_FASTQ_THREADS
    wrapper:
        SAMTOOLS_FASTQ_SEPARATE_WRAPPER
