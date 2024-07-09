import pandas as pd
from pathlib import Path

file_meta_df = pd.read_csv(snakemake.input[0], sep="\t").set_index(
    "file_id", drop=False, verify_integrity=True
)
for file_id in file_meta_df["file_id"]:
    Path.touch(f"{snakemake.output[0]}/{file_id}.bam")
