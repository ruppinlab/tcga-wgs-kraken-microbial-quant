import pandas as pd

use_cols = [
    "taxonomy_id",
    "kraken_assigned_reads",
    "added_reads",
    "new_est_reads",
]
idx_col_name = "taxonomy_id"

count_sum_df = pd.DataFrame()
for count_file in snakemake.input:
    count_df = pd.read_csv(
        count_file, sep="\t", header=0, index_col=idx_col_name, usecols=use_cols
    )[use_cols]
    count_sum_df = pd.concat([count_sum_df, count_df], axis=0).groupby(level=0).sum()

count_sum_df.index.name = idx_col_name
count_sum_df = count_sum_df.astype(int)
count_sum_df.to_csv(snakemake.output[0], sep="\t")
