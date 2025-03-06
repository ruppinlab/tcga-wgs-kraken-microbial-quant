# MAIT (Microbial Abundance In Tumors)

A Kraken2 + Bracken based pipeline for classifying and quantifying
microbial reads from GDC TCGA WGS data.

## Details

See README_details.md.

## Prerequistes

- This pipeline requires a token from the GDC that indicates access to the
  TCGA WGS read data.  Without such a token, the pipeline cannot proceed
  because the TCGA read data is controlled.

- The pipeline assumes a Unix or Linux-like environment.  It was tested
  on Red Hat Enterprise Linux 8.

- The pipeline analyzes > 10,000 samples and for practical use must be
  run on a cluster of computers.

- The complete pipeline uses around 4.5TB of disk space.

- The pipeline has only be tested on a cluster using the `slurm`
  scheduler (version 23.11.10).


## Installation

Install and set up a minimal conda enviroment using
[Miniforge3](https://github.com/conda-forge/miniforge#download)

The `conda` program will be used to install all required software.  Installation
was tested using conda 25.1.1, though any later version is expected to work.

Obtain the project source and create a conda environment with the tools
needed to run the project:

```bash
unzip MAITv1.0.0.zip
cd MAITv1.0.0
conda env create -f envs/tcga-wgs-kraken-microbial-quant.yaml
conda activate tcga-wgs-kraken-microbial-quant
```

This environment will, in particular contain Snakemake v8.16.

## Execution

Set the environment variable `GDC_TOKEN` to contain your GDC controlled
access authentication token.

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

## Dependencies

A Snakemake pipeline is divided into "rules" that run specific programs.  As
is typical for Snakemake pipelines, the pipeline rules will install the
appropriate versions of the software.   The dependent software is

  - samtools=1.20
  - hisat2=2.2.1
  - krakentools=1.2
  - kraken2=2.1.3
  - bracken=3.0
  - biobambam=2.0.183

Additionally, we customized a few files in the open-source Kraken 2 (v2.13),
Krakentools (v1.2), and Bracken (v3.0) projects.   The modified source, with
an explanation of what was changed, can be found in the `external/` directory
see `external/README_external.md`.  Note that MAIT has been licenced under
the GPL version 3.0.

## References

1. Lu et al. [Metagenome analysis using the Kraken software suite](
   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9725748/).
   Nat Protoc. 2022 Dec;17(12):2815-2839. doi: 10.1038/s41596-022-00738-y
2. Ge et al. [Comprehensive analysis of microbial content in whole-genome
   sequencing samples from The Cancer Genome Atlas project](
   https://doi.org/10.1101/2024.05.24.595788). bioRxiv 2024.05.24.595788
