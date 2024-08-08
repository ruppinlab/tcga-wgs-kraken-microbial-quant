from contextlib import redirect_stdout, redirect_stderr

with open(snakemake.log[0], "wt") as log_fh:
    with redirect_stdout(log_fh), redirect_stderr(log_fh):
        with open(snakemake.input.idmap, "rt") as id_ifh:
            with open(snakemake.output.idmap, "wt") as id_ofh:
                seqid2taxid_map = {}
                for line in id_ifh:
                    seqid, taxid = line.strip().replace(" ", "").split("\t", maxsplit=2)
                    seqid2taxid_map[seqid] = taxid
                    id_ofh.write(f"TAXID\tkraken:taxid|{taxid}|{seqid}\t{taxid}\n")
                with open(snakemake.input.fasta, "rt") as fa_ifh:
                    with open(snakemake.output.fasta, "wt") as fa_ofh:
                        for line in fa_ifh:
                            line = line.strip()
                            if line[0] == ">":
                                skip = False
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
                                    skip = True
                            if not skip:
                                fa_ofh.write(f"{line}\n")
