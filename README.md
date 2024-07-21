# tcga-wgs-kraken-microbial-quant

A Kraken-based pipeline for classifying and quantifying microbial
reads from GDC TCGA WGS data. Properly handles TCGA WGS merged
BAMs with mixed SE and PE reads and multiple read lengths when needed
by splitting to read group level FASTQs and processing data through
the pipeline at read group level before aggregating the results.
See [References](#references) for the general basis for this pipeline
and more information.

Pipeline short summary and full pipeline Snakemake rulegraph:

```
GDC TCGA WGS Unmapped Read BAMs ->
Biobambam2 Unmapped Read FASTQs (Split to Read Group Level when BAM has mixed PE/SE or mixed read lengths) ->
Bowtie2 Host Filtering (with T2T-CHM13v2.0) ->
Host Filtered BAMs ->
Biobambam2 Host Filtered FASTQs ->
Kraken2 Read Classification (with additional Protein Translated Search) ->
Bracken Read Quantification ->
Aggregate Read Group Level Counts ->
Count Matrix
```

![Snakemake rule graph](tcga-wgs-kraken-microbial-quant.svg)

## Prerequisites

The project was developed under GNU Linux and MacOS and assumes the
use of a Unix command line shell. Both Linux and MacOS provide a
command line shell by default. Other needed tools will be installed
by the instructions below.

## Installation

Install and set up
[Miniforge3](https://github.com/conda-forge/miniforge#miniforge3)

Obtain the project source and create a conda environment with the tools
needed to run the project:

```bash
git clone https://github.com/hermidalc/tcga-wgs-kraken-microbial-quant.git
cd tcga-wgs-kraken-microbial-quant
mamba env create -f envs/tcga-wgs-kraken-microbial-quant.yaml
mamba activate tcga-wgs-kraken-microbial-quant
```

Test that the installation is working by doing a dry run (if you don't
have a GDC token yet and wish to test your install do
`GDC_TOKEN='' snakemake --dry-run`). Below are the job statistics you
would see for all TCGA WGS primary tumors:

```
$ snakemake --dry-run

Building DAG of jobs...
Job stats:
job                          count
-------------------------  -------
all                              1
bowtie2_filtered_fastq_pe    13531
bowtie2_filtered_fastq_se      424
bowtie2_host_filter          13955
bowtie2_host_index               1
bracken_count_matrix             1
bracken_db_build                 7
bracken_merged_rg_counts        82
bracken_read_quant           13955
gdc_unmapped_bam             10838
gdc_unmapped_fastq_pe        13531
gdc_unmapped_fastq_se          424
host_genome_fasta                1
kraken2_db_build                 1
kraken2_db_library               5
kraken2_db_taxonomy              1
kraken2_read_classif         13955
total                        80713
```


## Execution

Given the compute intensive nature of this pipeline and the number of
jobs to execute it we highly recommend running the pipeline on an HPC
cluster.

Set your GDC controlled access authentication token in the environment
variable `GDC_TOKEN` or `GDC_TOKEN_FILE`, or the file `~/.gdc_token`
so the pipeline can get the token.

Run the workflow:

```bash
snakemake --use-conda --printshellcmds
```

## References

1. Lu et al. [Metagenome analysis using the Kraken software suite](
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9725748/).
Nat Protoc. 2022 Dec;17(12):2815-2839. doi: 10.1038/s41596-022-00738-y
2. Ge et al. [Comprehensive analysis of microbial content in whole-genome
sequencing samples from The Cancer Genome Atlas project](
    https://doi.org/10.1101/2024.05.24.595788). bioRxiv 2024.05.24.595788
