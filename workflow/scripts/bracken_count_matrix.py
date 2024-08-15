import pandas as pd

use_cols = ["taxonomy_id", "new_est_reads"]
idx_col_name = "taxonomy_id"

count_files = snakemake.input.get("counts") or snakemake.input
assert count_files is not None, "input: counts is a required parameter"
sample_names = snakemake.params.get("samples")
assert sample_names is not None, "params: samples is a required parameter"

count_matrix_df = pd.DataFrame()
for count_file, sample_name in zip(count_files, sample_names):
    count_df = pd.read_csv(
        count_file, sep="\t", header=0, index_col=idx_col_name, usecols=use_cols
    )[use_cols]
    count_df.columns = [sample_name]
    count_matrix_df = pd.concat(
        [count_matrix_df, count_df], axis=1, verify_integrity=True
    )

count_matrix_df.fillna(0, inplace=True)
count_matrix_df = count_matrix_df.astype(int)
count_matrix_df.index.name = idx_col_name
count_matrix_df.to_csv(snakemake.output[0], sep="\t")
