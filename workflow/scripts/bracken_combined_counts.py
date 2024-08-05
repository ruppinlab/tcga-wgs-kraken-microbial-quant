import pandas as pd

count_sum_df = pd.DataFrame()
for count_file in snakemake.input:
    count_df = pd.read_csv(
        count_file, sep="\t", header=0, index_col=0, usecols=list(range(6))
    )
    count_sum_df = pd.concat([count_sum_df, count_df], axis=0).groupby(level=0).sum()

count_sum_df.index.name = "name"
count_sum_df = count_sum_df.astype(int)
count_sum_df.to_csv(snakemake.output[0], sep="\t")
