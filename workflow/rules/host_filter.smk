from os.path import dirname


rule bowtie2_index:
    input:
        ref=REF_FASTA_FILE,
    params:
        extra=f"--seed {config['random_seed']}",
    output:
        BOWTIE2_INDEX_FILES,
    log:
        BOWTIE2_INDEX_LOG,
    threads: BOWTIE2_BUILD_THREADS
    wrapper:
        BOWTIE2_BUILD_WRAPPER


rule bowtie2_align_pe:
    input:
        sample=[GDC_UNMAPPED_FASTQ_R1_FILE, GDC_UNMAPPED_FASTQ_R2_FILE],
        idx=BOWTIE2_INDEX_FILES,
    params:
        extra=f"--seed {config['random_seed']}",
    output:
        temp(BOWTIE2_SAM_PE_FILE),
        unaligned=BOWTIE2_FILTERED_SAM_PE_FILE,
    log:
        BOWTIE2_ALIGN_LOG,
    threads: BOWTIE2_ALIGN_THREADS
    wrapper:
        BOWTIE2_ALIGN_WRAPPER


rule bowtie2_align_se:
    input:
        sample=GDC_UNMAPPED_FASTQ_SE_FILE,
        idx=BOWTIE2_INDEX_FILES,
    params:
        extra=f"--seed {config['random_seed']}",
    output:
        temp(BOWTIE2_SAM_SE_FILE),
        unaligned=BOWTIE2_FILTERED_SAM_SE_FILE,
    log:
        BOWTIE2_ALIGN_LOG,
    threads: BOWTIE2_ALIGN_THREADS
    wrapper:
        BOWTIE2_ALIGN_WRAPPER


rule bowtie2_sorted_sam_pe:
    input:
        BOWTIE2_FILTERED_SAM_PE_FILE,
    params:
        extra=config["samtools"]["sort"]["extra"],
    output:
        BOWTIE2_SORTED_FILTERED_SAM_PE_FILE,
    log:
        BOWTIE2_SORTED_FILTERED_SAM_LOG,
    threads: SAMTOOLS_SORT_THREADS
    wrapper:
        SAMTOOLS_SORT_WRAPPER


rule bowtie2_sorted_sam_se:
    input:
        BOWTIE2_FILTERED_SAM_SE_FILE,
    params:
        extra=config["samtools"]["sort"]["extra"],
    output:
        BOWTIE2_SORTED_FILTERED_SAM_SE_FILE,
    log:
        BOWTIE2_SORTED_FILTERED_SAM_LOG,
    threads: SAMTOOLS_SORT_THREADS
    wrapper:
        SAMTOOLS_SORT_WRAPPER


rule bowtie2_filtered_fastq_pe:
    input:
        BOWTIE2_SORTED_FILTERED_SAM_PE_FILE,
    params:
        extra=config["samtools"]["fastq"]["extra"],
    output:
        BOWTIE2_FILTERED_FASTQ_R1_FILE,
        BOWTIE2_FILTERED_FASTQ_R2_FILE,
    log:
        BOWTIE2_FILTERED_FASTQ_LOG,
    threads: SAMTOOLS_FASTQ_THREADS
    wrapper:
        SAMTOOLS_FASTQ_SEPARATE_WRAPPER


rule bowtie2_filtered_fastq_se:
    input:
        BOWTIE2_SORTED_FILTERED_SAM_SE_FILE,
    params:
        extra=config["samtools"]["fastq"]["extra"],
    output:
        BOWTIE2_FILTERED_FASTQ_SE_FILE,
    log:
        BOWTIE2_FILTERED_FASTQ_LOG,
    threads: SAMTOOLS_FASTQ_THREADS
    wrapper:
        SAMTOOLS_FASTQ_SEPARATE_WRAPPER
