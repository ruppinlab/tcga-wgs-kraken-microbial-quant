import pandas as pd

cols = [
    "length",
    "reads",
    "pct_reads",
    "cum_reads",
    "cum_pct_reads",
    "bases",
    "pct_bases",
    "cum_bases",
    "cum_pct_bases",
]
print(f"Getting maximum read length from {snakemake.input[0]}", flush=True)
df = pd.read_csv(snakemake.input[0], sep="\t", comment="#", names=cols)
max_length = df["length"].max().astype(str)
with open(snakemake.output[0], "w") as fh:
    fh.write(max_length)
with open(snakemake.log[0], "w") as fh:
    fh.write(f"Wrote max length {max_length} file")
