# tcga-wgs-kraken-microbial-quant

A Kraken-based pipeline for classifying and quantifying microbial
reads from GDC TCGA WGS data. Supports Kraken2 and KrakenUniq for read
classification. Includes the option to do a second pass Kraken2 protein
translated search of the unclassified reads from the Kraken2 first pass
and combining the report results before feeding into Bracken. Properly
handles TCGA WGS merged BAMs with mixed SE and PE reads and multiple
read lengths when needed by splitting them into read group level FASTQs
and processing data through the pipeline at read group level before
aggregating the Bracken count results back to GDC BAM level.

See [References](#references) for the general basis for this pipeline
and more information.

A high-level pipeline summary:

```
GDC TCGA WGS Unmapped Read BAMs ->
Biobambam2 Unmapped Read FASTQs (Split to Read Group Level when BAM has mixed PE/SE or read lengths) ->
Bowtie2 Host Filtering (with T2T-CHM13v2.0) ->
Biobambam2 Host Filtered FASTQs ->
Kraken2 Nucleotide Read Classification ->
Kraken2 Translated Protein Search Read Classification of Unclassified Reads ->
KrakenTools Combine Nucleotide and Protein Reports ->
Bracken Read Quantification ->
Aggregate Read Group Level Counts ->
Count Matrix
```

Full Snakemake pipeline
rulegraph:

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
would see for an analysis of all TCGA WGS primary tumors using Kraken2
with the optional addtional step of Kraken2 protein translated search of
unclassified reads:

```
$ snakemake --dry-run

Building DAG of jobs...
Job stats:
job                             count
----------------------------  -------
all                                 1
bowtie2_filtered_fastq_pe       13531
bowtie2_filtered_fastq_se         424
bowtie2_host_filter             13955
bowtie2_host_index                  1
bracken_count_matrix                1
bracken_db_build                    7
bracken_merged_rg_counts           82
bracken_read_quant              13955
gdc_unmapped_bam                10838
gdc_unmapped_fastq_pe           13531
gdc_unmapped_fastq_se             424
host_genome_fasta                   1
kraken2_combined_report         13955
kraken2_db_build                    2
kraken2_db_taxonomy                 2
kraken2_nucl_library                5
kraken2_nucl_read_classif_pe    13531
kraken2_nucl_read_classif_se      424
kraken2_prot_library                4
kraken2_prot_read_classif_pe    13531
kraken2_prot_read_classif_se      424
total                          108629
```


## Execution

Given the compute intensive nature of this pipeline and the large
number of jobs to execute it we highly recommend running the it on
an HPC cluster.

Set your GDC controlled access authentication token in the environment
variable `GDC_TOKEN` or `GDC_TOKEN_FILE`, or the file `~/.gdc_token`
so the pipeline can get the token.

Run the workflow:

```bash
snakemake --use-conda --printshellcmds
```

The pipeline is configured to not require much storage, as intermediate
files are flagged as temporary and deleted when they are no longer
needed as the pipeline is running, except for Bracken outputs. If you
would like to keep intermediate files for other uses, specifiy the
`--notemp` snakemake option in the execution command above.

## References

1. Lu et al. [Metagenome analysis using the Kraken software suite](
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9725748/).
Nat Protoc. 2022 Dec;17(12):2815-2839. doi: 10.1038/s41596-022-00738-y
2. Ge et al. [Comprehensive analysis of microbial content in whole-genome
sequencing samples from The Cancer Genome Atlas project](
    https://doi.org/10.1101/2024.05.24.595788). bioRxiv 2024.05.24.595788
