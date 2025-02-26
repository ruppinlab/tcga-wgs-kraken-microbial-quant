# External Tools

This directory contains external open-source tools that were modified
and are used in this project.  We modified parts of Bracken, kraken2
and KrakenTools.

Each subdirectory contains an archive of the unmodified source of the
named open-source project, and a "custom" directory containing an
unmodified LICENCE and only the source files that were modified.  For
example, the kraken2 has the following structure:

    kraken2
    ├── 2.1.3-custom
    │   ├── LICENSE
    │   └── scripts
    │       └── k2
    └── 2.1.3.tar.gz

As shown by this tree, in kraken2, only the scripts/k2 file was 
modified.

We provide the external source as an archive to indicate that it has
not been modified.  We include an unmodified license file with the
modified source in an abundance of caution; the same license file is
in the archive.

We obtained archives of the original projects by archiving a working
copy of the project checked out at the specified tag.

https://github.com/DerrickWood/kraken2.git        tag: v2.1.3
https://github.com/jenniferlu717/Bracken.git      tag: v3.0
https://github.com/jenniferlu717/KrakenTools.git  tag: v2.0

We make the following changes.

Kraken2
=======
We customized ‘k2’ program from the open-source MIT licensed Kraken
2.1.3 distribution.  The customizations make the script compatible
with the unique environment of the NIH Biowulf cluster, though the
code changes do not restrict use of the script to Biowulf.

Bracken
=======
We customized the open-source GPL licensed Bracken 3.0 software to
make it optional a lower bound on the unique k-mer count when deciding
to whether calls matched to a specific taxonomy level are spurious.

KrakenTools
===========
We modified `combine_kreports.py` to improve handling of empty trees
which are common in our pipeline.
