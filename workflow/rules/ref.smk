rule ref_genome_seq:
    params:
        REF_FASTA_URL,
    output:
        REF_FASTA_FILE,
    message:
        "{params}"
    retries: config["download"]["retries"]
    script:
        "../scripts/url_file.py"
