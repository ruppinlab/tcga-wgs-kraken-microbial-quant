import re
from contextlib import redirect_stdout, redirect_stderr

import pandas as pd

with open(snakemake.log[0], "wt") as log_fh:
    with redirect_stdout(log_fh), redirect_stderr(log_fh):
        with open(snakemake.output.idmap, "wt") as id_ofh:
            with open(snakemake.output.fasta, "wt") as fa_ofh:
                with open(snakemake.input.fasta, "rt") as fa_ifh:
                    skipped_organisms = []
                    meta_df = pd.read_csv(snakemake.input.meta, sep="\t")
                    for line in fa_ifh:
                        line = line.strip()
                        if line[0] == ">":
                            skip = False
                            organism = re.findall(
                                r"\s*\|\s*organism=(.+?)\s*\|\s*", line, re.IGNORECASE
                            )
                            organism = organism[0].replace("_", " ")
                            if organism not in meta_df["Organism"].values:
                                skip = True
                                if organism not in skipped_organisms:
                                    print(f"{organism} metadata not found, skipping")
                                    skipped_organisms.append(organism)
                            line_parts = re.split(r"\s*\|\s*", line)
                            seqid = line_parts[0].lstrip(">")
                            taxid = meta_df.loc[
                                meta_df["Organism"] == organism, "Species NCBI taxon ID"
                            ].squeeze()
                            id_ofh.write(
                                f"TAXID\tkraken:taxid|{taxid}|{seqid}\t{taxid}\n"
                            )
                            line = f">kraken:taxid|{taxid}|{seqid}"
                        if not skip:
                            fa_ofh.write(f"{line}\n")
