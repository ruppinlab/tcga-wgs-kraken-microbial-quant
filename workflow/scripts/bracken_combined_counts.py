import pandas as pd

idx_cols = ["name", "taxonomy_id", "taxonomy_lvl"]

count_dfs = []
np_pair = snakemake.params.get("np_pair", False)

if np_pair:
    if len(snakemake.input) != 2:
        raise ValueError(
            "Input must be a nucl and prot pair of count files when --np-pair is set"
        )
    nucl_df = pd.read_csv(snakemake.input[0], sep="\t", header=0, index_col=idx_cols)
    prot_df = pd.read_csv(snakemake.input[1], sep="\t", header=0, index_col=idx_cols)
    if not nucl_df.empty:
        nucl_df.drop(columns=["fraction_total_reads"], inplace=True)
        count_dfs.append(nucl_df)
        prot_df = prot_df.loc[prot_df.index.isin(nucl_df.index)]
        if not prot_df.empty:
            prot_df.drop(columns=["fraction_total_reads"], inplace=True)
            count_dfs.append(prot_df)
else:
    for count_file in snakemake.input:
        count_df = pd.read_csv(count_file, sep="\t", header=0, index_col=idx_cols)
        if not count_df.empty:
            count_df.drop(columns=["fraction_total_reads"], inplace=True)
            count_dfs.append(count_df)

if count_dfs:
    all_count_df = pd.concat(count_dfs, axis=0)
    count_sum_df = all_count_df.groupby(idx_cols).sum()
    count_sum_df["fraction_total_reads"] = round(
        count_sum_df["new_est_reads"] / count_sum_df["new_est_reads"].sum(), 5
    )
else:
    count_sum_df = pd.DataFrame(
        columns=[
            "name",
            "taxonomy_id",
            "taxonomy_lvl",
            "kraken_assigned_reads",
            "added_reads",
            "new_est_reads",
            "fraction_total_reads",
        ]
    )

count_sum_df.to_csv(snakemake.output[0], sep="\t")
