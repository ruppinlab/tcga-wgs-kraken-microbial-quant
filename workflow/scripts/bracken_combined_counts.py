import pandas as pd

idx_cols = ["name", "taxonomy_id", "taxonomy_lvl"]

count_dfs = []
for count_file in snakemake.input:
    count_df = pd.read_csv(count_file, sep="\t", header=0, index_col=idx_cols)
    count_df.drop(columns=["fraction_total_reads"], inplace=True)
    if not count_df.empty:
        count_dfs.append(count_df)

all_count_df = pd.concat(count_dfs, axis=0)
count_sum_df = all_count_df.groupby(idx_cols).sum()
count_sum_df["fraction_total_reads"] = round(
    count_sum_df["new_est_reads"] / count_sum_df["new_est_reads"].sum(), 5
)
count_sum_df.to_csv(snakemake.output[0], sep="\t")
