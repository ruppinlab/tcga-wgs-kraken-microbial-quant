import pandas as pd

sum_count_df = pd.DataFrame()

for count_file in snakemake.input:
    count_df = pd.read_csv(
        count_file, sep="\t", header=0, index_col=0, usecols=list(range(6))
    )
    sum_count_df = pd.concat([sum_count_df, count_df], axis=0).groupby(level=0).sum()

sum_count_df.index.name = "name"
sum_count_df = sum_count_df.astype(int)
sum_count_df.sort_index(inplace=True)
sum_count_df.to_csv(snakemake.output[0], sep="\t")
