# MAIT (Microbial Abundance In Tumors)

A Kraken2 + Bracken based pipeline for classifying and quantifying
microbial reads from GDC TCGA WGS data.

## Details

See README_details.md.

## Prerequistes

- This pipeline requires a GDC token indicating access to the TCGA WGS
  read data.  Without such a token, the pipeline cannot proceed
  because the TCGA read data is controlled.

- The pipeline analyzes > 10,000 samples and for practical use must be
  run on a cluster of computers.

- The complete pipeline uses around 4.5TB of disk space.

- The pipeline has only be tested on a cluster using the `slurm`
  scheduler.

## Installation

Install and set up
[Miniforge3](https://github.com/conda-forge/miniforge#download)

Obtain the project source and create a conda environment with the tools
needed to run the project:

```bash
unzip MAITv1.0.0.zip
cd MAITv1.0.0
mamba env create -f envs/tcga-wgs-kraken-microbial-quant.yaml
mamba activate tcga-wgs-kraken-microbial-quant
```

## Execution

Set your GDC controlled access authentication token in the environment
variable `GDC_TOKEN`

Run the workflow on a cluster:

```bash
./scripts/run_snakemake_slurm.sh \
--workflow-profile workflow/profiles/biowulf \
--sbatch-opts="--time=3-00:00:00 --cpus-per-task=28 --mem=10248"
```

We've provided a SLURM cluster configuration for the NIH HPC cluster,
though it is straightforward to create a Snakemake v8+ cluster config for
your particular needs.

In principle, the workflow can be adapted to any system for which
`snakemake` can submit jobs, but we have only tested in on slurm.

## Output

The raw counts for each sample by each genus will be found in

`results/bracken/matrix/tcga_wgs_primary_tumor_genus_count_matrix.tsv`

## Post-processing

The directory `manuscript_code_and_data` contains code for
post-processing the `tcga_wgs_primary_tumor_genus_count_matrix.tsv`
and producing tables of raw counts and CPM data.  See README.md file
in the subdirectory for instructions.

## References

1. Lu et al. [Metagenome analysis using the Kraken software suite](
   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9725748/).
   Nat Protoc. 2022 Dec;17(12):2815-2839. doi: 10.1038/s41596-022-00738-y
2. Ge et al. [Comprehensive analysis of microbial content in whole-genome
   sequencing samples from The Cancer Genome Atlas project](
   https://doi.org/10.1101/2024.05.24.595788). bioRxiv 2024.05.24.595788
