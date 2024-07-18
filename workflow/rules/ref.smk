rule host_genome_fasta:
    params:
        HOST_REF_FASTA_URL,
    output:
        HOST_REF_FASTA_FILE,
    message:
        "{params}"
    retries: config["download"]["retries"]
    script:
        "../scripts/url_file.py"
