import pandas as pd

count_files = snakemake.input.get("counts") or snakemake.input
assert count_files is not None, "input: counts is a required parameter"
sample_names = snakemake.params.get("samples")
assert sample_names is not None, "params: samples is a required parameter"

count_matrix_df = pd.DataFrame()
for count_file, sample_name in zip(count_files, sample_names):
    count_df = pd.read_csv(
        count_file, sep="\t", header=None, index_col=0, usecols=[0, 3]
    )
    count_df.columns = [sample_name]
    count_matrix_df = pd.concat(
        [count_matrix_df, count_df], axis=1, verify_integrity=True
    )

count_matrix_df = count_matrix_df.astype(int)
count_matrix_df.fillna(0)
count_matrix_df.index.name = "Taxonomic Name"
count_matrix_df.to_csv(snakemake.output[0], sep="\t")
