rule eupathdb_fasta_archive:
    params:
        EUPATHDB_FASTA_TGZ_URL,
    output:
        temp(EUPATHDB_FASTA_TGZ_FILE),
    log:
        EUPATHDB_FASTA_TGZ_LOG,
    message:
        "{params}"
    retries: config["download"]["retries"]
    conda:
        "../envs/wget.yaml"
    shell:
        "wget -nv -O {output} {params} > {log} 2>&1"


rule eupathdb_fastas:
    input:
        EUPATHDB_FASTA_TGZ_FILE,
    params:
        EUPATHDB_RESOURCES_DIR,
    output:
        temp(directory(EUPATHDB_FASTA_DIR)),
    log:
        EUPATHDB_FASTA_LOG,
    shell:
        "tar -xvzf {input} -C {params} > {log} 2>&1"


rule eupathdb_merged_fasta:
    input:
        EUPATHDB_FASTA_DIR,
    output:
        EUPATHDB_MERGED_FASTA_FILE,
    log:
        EUPATHDB_MERGED_FASTA_LOG,
    shell:
        "find {input} -type f -name '*.fna' -exec cat {{}} + 1> {output} 2> {log}"


rule eupathdb_seqid2taxid_map:
    params:
        EUPATHDB_SEQID2TAXID_MAP_URL,
    output:
        temp(EUPATHDB_SEQID2TAXID_MAP_FILE),
    log:
        EUPATHDB_SEQID2TAXID_MAP_LOG,
    message:
        "{params}"
    retries: config["download"]["retries"]
    conda:
        "../envs/wget.yaml"
    shell:
        "wget -nv -O {output} {params} > {log} 2>&1"


rule kraken2_eupathdb_lib_files:
    input:
        idmap=EUPATHDB_SEQID2TAXID_MAP_FILE,
        fasta=EUPATHDB_MERGED_FASTA_FILE,
    output:
        idmap=KRAKEN2_EUPATHDB_LIB_IDMAP_FILE,
        fasta=KRAKEN2_EUPATHDB_LIB_FASTA_FILE,
    log:
        KRAKEN2_EUPATHDB_LIB_FASTA_LOG,
    run:
        from contextlib import redirect_stdout, redirect_stderr

        with open(log[0], "wt") as log_fh:
            with redirect_stdout(log_fh), redirect_stderr(log_fh):
                seqid2taxid_map = {}
                with open(output.idmap, "wt") as id_ofh:
                    with open(input.idmap, "rt") as id_ifh:
                        for line in id_ifh:
                            seqid, taxid = (
                                line.strip().replace(" ", "").split("\t", maxsplit=2)
                            )
                            seqid2taxid_map[seqid] = taxid
                            id_ofh.write(
                                f"TAXID\tkraken:taxid|{taxid}|{seqid}\t{taxid}\n"
                            )
                    with open(output.fasta, "wt") as fa_ofh:
                        with open(input.fasta, "rt") as fa_ifh:
                            for line in fa_ifh:
                                line = line.strip()
                                if line[0] == ">":
                                    seqid = (
                                        line.replace(" ", "")
                                        .replace(">", "")
                                        .replace("|", "")
                                    )
                                    if seqid in seqid2taxid_map:
                                        taxid = seqid2taxid_map[seqid]
                                        line = f">kraken:taxid|{taxid}|{seqid}"
                                    elif seqid.startswith("LMARLEM2494"):
                                        # fix for missing Leishmania martiniquensis LEM2494
                                        taxid = 1580590
                                        line = f">kraken:taxid|{taxid}|{seqid}"
                                        id_ofh.write(
                                            f"TAXID\tkraken:taxid|{taxid}|{seqid}\t{taxid}\n"
                                        )
                                    else:
                                        log_fh.write(f"No taxid: {seqid}\n")
                                fa_ofh.write(f"{line}\n")
