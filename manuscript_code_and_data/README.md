# Code and Data

Code to post-process the output of MAIT (Microbial Abundance In Tumors)
and reproduce the tables and figures published in 

    Gertz et al. A revised estimation of microbial abundances in the Tumor
    Cancer Genome Atlas.

Because MAIT is a long-running pipeline, but scripts for
post-processing might be independently useful, we provide the output
file produced in our run of MAIT and used in our manuscript as

     tcga_wgs_primary_tumor_genus_count_matrix.tsv

To run the post-processing using a different raw count matrix, simply
replace this file with the output of MAIT.

The code is provided as R markdown (.Rmd) files that are most
conveniently run using Rstudio
(https://posit.co/download/rstudio-desktop/).  The resources of most
modern laptops or desktops (circa 2025) should easily be sufficient to
run these scripts.  We provide the code as R markdown documents
because each R markdown script displays intermediate results
graphically, allows the user to inspect intermediate results, and
contains text that documents its operation. 

We remark that we also provide, within this directory, the output of
all scripts.  Therefore users who wish to modify nothing may just take
the output files as given without running the scripts.

## Prerequiste libraries

Install the following libraries using the `install.packages` function.
```
c('assert', 'dendextend', 'janitor', 'pheatmap', 
'plotrix', 'readxl', 'tidyverse', 'vegan')
```
## Post-processing reads

Run the R markdown file `process_pipeline_output.Rmd`.  The pipeline
output used in the paper is already provided in this directory in the file

    tcga_wgs_primary_tumor_genus_count_matrix.tsv

Replace this file if you made any changes to the pipeline that may
have changed the output.

The `process_pipeline_output.Rmd` script also uses as input the files

    tcga_wgs_primary_tumor_file_meta.tsv
    total_reads.tsv

provided in this directory.  Its output are the files

- `S1_filtered_cpm.tsv`
- `S4_filtered_counts.tsv`
- `S5_raw_counts.tsv`
- `S6_raw_cpm.tsv`

Which represent Tables S1, S4, S4 and S6 of the manuscript.

## Generate figures

Run `generate_figures.Rmd`.  This script reads `S1_filtered_cpm.tsv`
and `tcga_wgs_primary_tumor_file_meta.tsv`.  It produces Figures 1B,
2A-D.  Figure files should be visually the same as those provided
in this archive. 

Some differences in the bitwise content (changes in metadata) are
expected each time the figures are generated.  It is not an error
that the output of the script are not bitwise identical to the files
in the archive.

The final figures presented in the manuscript were subject to
noncontroversial manual cosmetic changes.

## Search databases for human-associated microbes

The script `identify_human_host.Rmd` contains code to search several
microbial databases for microbes known to be seen in human.
Unfortunately, the databases themselves are not redistributable and
may be behind paywalls.

There are instructions at the top of `identify_human_host.Rmd` that
tell you how to retrieve the data and what filenames and formats are
expected by the script.

The result of the script is `raw_human_host.tsv` which lists the
genera and the evidence, if any, from the database search that each
genus has been observed in human.

Hand curation by literature search is beyond the scope of this code;
see the manuscript and Tables S2 and S3 for the result of hand
curation.
