import pandas as pd

use_cols = [
    "taxonomy_id",
    "kraken_assigned_reads",
    "added_reads",
    "new_est_reads",
]

all_count_df = pd.DataFrame()
for count_file in snakemake.input:
    count_df = pd.read_csv(
        count_file, sep="\t", header=0, index_col=use_cols[0], usecols=use_cols
    )[use_cols[1:]]
    if not count_df.empty:
        all_count_df = pd.concat([all_count_df, count_df], axis=0)

count_sum_df = all_count_df.groupby(level=0).sum()
count_sum_df.index.name = use_cols[0]
count_sum_df = count_sum_df.astype(int)
count_sum_df.to_csv(snakemake.output[0], sep="\t")
