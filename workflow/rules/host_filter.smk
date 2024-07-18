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


rule bowtie2_align:
    input:
        sample=lambda wc: (
            [GDC_UNMAPPED_FASTQ_R1_FILE, GDC_UNMAPPED_FASTQ_R2_FILE]
            if wc.etype == "pe"
            else GDC_UNMAPPED_FASTQ_SE_FILE
        ),
        idx=BOWTIE2_INDEX_FILES,
    params:
        extra=f"--seed {config['random_seed']}",
    output:
        temp(BOWTIE2_SAM_FILE),
        unaligned=BOWTIE2_FILTERED_SAM_FILE,
    log:
        BOWTIE2_ALIGN_LOG,
    threads: BOWTIE2_ALIGN_THREADS
    wrapper:
        BOWTIE2_ALIGN_WRAPPER


rule bowtie2_filtered_fastq_pe:
    input:
        BOWTIE2_FILTERED_SAM_FILE,
    params:
        paired_end=True,
        extra=(
            f"{config['biobambam2']['bamtofastq']['extra']['common']} "
            f"{config['biobambam2']['bamtofastq']['extra']['paired_end']}"
        ),
    output:
        BOWTIE2_FILTERED_FASTQ_R1_FILE,
        BOWTIE2_FILTERED_FASTQ_R2_FILE,
    log:
        BOWTIE2_FILTERED_FASTQ_LOG,
    wrapper:
        BIOBAMBAM2_BAMTOFASTQ_WRAPPER


rule bowtie2_filtered_fastq_se:
    input:
        BOWTIE2_FILTERED_SAM_FILE,
    params:
        paired_end=False,
        extra=config["biobambam2"]["bamtofastq"]["extra"]["common"],
    output:
        BOWTIE2_FILTERED_FASTQ_SE_FILE,
    log:
        BOWTIE2_FILTERED_FASTQ_LOG,
    wrapper:
        BIOBAMBAM2_BAMTOFASTQ_WRAPPER
