rule host_genome_fasta:
    params:
        HOST_REF_FASTA_URL,
    output:
        HOST_REF_FASTA_FILE,
    message:
        "{params}"
    retries: config["download"]["retries"]
    run:
        from urllib.request import urlcleanup, urlretrieve

        urlretrieve(snakemake.params[0], filename=snakemake.output[0])
        urlcleanup()


rule bowtie2_host_index:
    input:
        ref=HOST_REF_FASTA_FILE,
    params:
        extra=f"--seed {config['random_seed']}",
    output:
        BOWTIE2_HOST_INDEX_FILES,
    log:
        BOWTIE2_HOST_INDEX_LOG,
    threads: BOWTIE2_BUILD_THREADS
    wrapper:
        BOWTIE2_BUILD_WRAPPER


rule bowtie2_host_filter:
    input:
        sample=lambda wc: (
            [GDC_UNMAPPED_FASTQ_R1_FILE, GDC_UNMAPPED_FASTQ_R2_FILE]
            if wc.etype == "pe"
            else GDC_UNMAPPED_FASTQ_SE_FILE
        ),
        idx=BOWTIE2_HOST_INDEX_FILES,
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
        BOWTIE2_FILTERED_SAM_PE_FILE,
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
        BOWTIE2_FILTERED_FASTQ_PE_LOG,
    wrapper:
        BIOBAMBAM2_BAMTOFASTQ_WRAPPER


rule bowtie2_filtered_fastq_se:
    input:
        BOWTIE2_FILTERED_SAM_SE_FILE,
    params:
        paired_end=False,
        extra=config["biobambam2"]["bamtofastq"]["extra"]["common"],
    output:
        BOWTIE2_FILTERED_FASTQ_SE_FILE,
    log:
        BOWTIE2_FILTERED_FASTQ_SE_LOG,
    wrapper:
        BIOBAMBAM2_BAMTOFASTQ_WRAPPER
