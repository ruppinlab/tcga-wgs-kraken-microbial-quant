rule host_genome_fasta:
    params:
        url=HOST_REF_FASTA_URL,
    output:
        HOST_REF_FASTA_FILE,
    log:
        HOST_REF_FASTA_LOG,
    message:
        "{params.url}"
    retries: config["download"]["retries"]
    conda:
        "../envs/wget.yaml"
    shell:
        "wget -nv -O - '{params.url}' | gunzip -c 1> {output} 2> {log}"


rule host_genome_index:
    input:
        fasta=HOST_REF_FASTA_FILE,
    params:
        prefix=HOST_GENOME_INDEX_PREFIX,
        extra=f"{config[HOST_FILTER_MODE]['build']['extra']} --seed {config['random_seed']}",
    output:
        directory(HOST_GENOME_INDEX_DIR),
    log:
        HOST_GENOME_INDEX_LOG,
    threads: HOST_BUILD_THREADS
    wrapper:
        HOST_BUILD_WRAPPER


rule host_filtered_fastq_pe:
    input:
        reads=[GDC_UNMAPPED_FASTQ_R1_FILE, GDC_UNMAPPED_FASTQ_R2_FILE],
        dir=HOST_GENOME_INDEX_DIR,
    params:
        idx=HOST_GENOME_INDEX_PREFIX,
        extra=f"{config[HOST_FILTER_MODE]['align']['extra']} --seed {config['random_seed']}",
    output:
        temp(HOST_BAM_PE_FILE),
        unconcordant=[
            temp(HOST_FILTERED_FASTQ_R1_FILE),
            temp(HOST_FILTERED_FASTQ_R2_FILE),
        ],
    log:
        HOST_FILTERED_FASTQ_PE_LOG,
    threads: HOST_ALIGN_THREADS
    wrapper:
        HOST_ALIGN_WRAPPER


rule host_filtered_fastq_se:
    input:
        reads=[GDC_UNMAPPED_FASTQ_SE_FILE],
        dir=HOST_GENOME_INDEX_DIR,
    params:
        idx=HOST_GENOME_INDEX_PREFIX,
        extra=f"{config[HOST_FILTER_MODE]['align']['extra']} --seed {config['random_seed']}",
    output:
        temp(HOST_BAM_SE_FILE),
        unaligned=temp(HOST_FILTERED_FASTQ_SE_FILE),
    log:
        HOST_FILTERED_FASTQ_SE_LOG,
    threads: HOST_ALIGN_THREADS
    wrapper:
        HOST_ALIGN_WRAPPER
