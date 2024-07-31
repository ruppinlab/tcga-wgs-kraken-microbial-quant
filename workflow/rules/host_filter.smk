rule host_genome_fasta:
    params:
        HOST_REF_FASTA_URL,
    output:
        HOST_REF_FASTA_FILE,
    log:
        HOST_REF_FASTA_LOG,
    message:
        "{params}"
    retries: config["download"]["retries"]
    conda:
        "../envs/wget.yaml"
    shell:
        "wget -nv -O {output} {params} > {log} 2>&1"


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


rule bowtie2_filtered_fastq_pe:
    input:
        sample=[GDC_UNMAPPED_FASTQ_R1_FILE, GDC_UNMAPPED_FASTQ_R2_FILE],
        idx=BOWTIE2_HOST_INDEX_FILES,
    params:
        extra=f"--seed {config['random_seed']}",
    output:
        temp(BOWTIE2_SAM_PE_FILE),
        un_conc=[
            temp(BOWTIE2_FILTERED_FASTQ_R1_FILE),
            temp(BOWTIE2_FILTERED_FASTQ_R2_FILE),
        ],
    log:
        BOWTIE2_FILTERED_FASTQ_PE_LOG,
    threads: BOWTIE2_ALIGN_THREADS
    wrapper:
        BOWTIE2_ALIGN_WRAPPER


rule bowtie2_filtered_fastq_se:
    input:
        sample=[GDC_UNMAPPED_FASTQ_SE_FILE],
        idx=BOWTIE2_HOST_INDEX_FILES,
    params:
        extra=f"--seed {config['random_seed']}",
    output:
        temp(BOWTIE2_SAM_SE_FILE),
        unaligned=temp(BOWTIE2_FILTERED_FASTQ_SE_FILE),
    log:
        BOWTIE2_FILTERED_FASTQ_SE_LOG,
    threads: BOWTIE2_ALIGN_THREADS
    wrapper:
        BOWTIE2_ALIGN_WRAPPER
