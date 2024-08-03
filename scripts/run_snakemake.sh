#!/bin/bash

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
while [[ -v CONDA_DEFAULT_ENV ]]; do
    conda deactivate
done
conda activate tcga-wgs-kraken-microbial-quant

SCRIPT_DIR=$(dirname $(realpath -s $0))
PROJECT_DIR=$(realpath $SCRIPT_DIR/../)

chdir $PROJECT_DIR

SNAKEMAKE_CMD="snakemake $@"
echo $SNAKEMAKE_CMD
$SNAKEMAKE_CMD
