__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

import re

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

index = snakemake.input.get("index")
assert index is not None, "input: index is a required input parameter"

n = len(snakemake.input.reads)
assert (
    n == 1 or n == 2
), "input: reads must have 1 (single-end) or 2 (paired-end) elements."

reads = ""
if n == 1:
    if get_format(snakemake.input.reads[0]) in ("bam", "sam"):
        reads = f"-b {snakemake.input.reads}"
    else:
        if snakemake.params.get("interleaved", False):
            reads = f"--interleaved {snakemake.input.reads}"
        else:
            reads = f"-U {snakemake.input.reads}"
else:
    reads = "-1 {} -2 {}".format(*snakemake.input.reads)

extra = snakemake.params.get("extra", "")
if all(get_format(r) in ("fastq", "fq") for r in snakemake.input.reads):
    extra += " -q "
elif all(get_format(r) == "tab5" for r in snakemake.input.reads):
    extra += " --tab5 "
elif all(get_format(r) == "tab6" for r in snakemake.input.reads):
    extra += " --tab6 "
elif all(get_format(r) in ("fa", "mfa", "fasta") for r in snakemake.input.reads):
    extra += " -f "

unaligned = snakemake.output.get("unaligned")
if unaligned:
    un = (
        "un-gz"
        if unaligned.endswith((".gz", ".GZ"))
        else "un-bz2" if unaligned.endswith((".bz2", ".BZ2")) else "un"
    )
    extra += f" --{un} {unaligned} "

unpaired = snakemake.output.get("unpaired")
if unpaired:
    al = (
        "al-gz"
        if unpaired.endswith((".gz", ".GZ"))
        else "al-bz2" if unpaired.endswith((".bz2", ".BZ2")) else "al"
    )
    extra += f" --{al} {unpaired} "

unconcordant = snakemake.output.get("unconcordant")
if unconcordant:
    assert (
        isinstance(unconcordant, (list, tuple, set)) and len(unconcordant) == 2
    ), "input: unconcordant must be a list of two files"
    un_conc_groups_1 = re.split(
        r"(?:1|2)\.(fq|fastq)(?:\.gz)?$", unconcordant[0], flags=re.IGNORECASE
    )
    un_conc_groups_2 = re.split(
        r"(?:1|2)\.(fq|fastq)(?:\.gz)?$", unconcordant[1], flags=re.IGNORECASE
    )
    assert set(un_conc_groups_1) == set(
        un_conc_groups_2
    ), "input: unconcordant file name pair is invalid"
    un_conc_opt = (
        "un-conc-gz"
        if un_conc_groups_1[2].lower() == ".gz"
        else "un-conc-bz2" if un_conc_groups_1[2].lower() == ".bz2" else "un-conc"
    )
    un_conc_pat = f"{un_conc_groups_1[0]}%.{un_conc_groups_1[1]}"
    extra += f" --{un_conc_opt} {un_conc_pat} "

concordant = snakemake.output.get("concordant")
if concordant:
    assert (
        isinstance(concordant, (list, tuple, set)) and len(concordant) == 2
    ), "input: concordant must be a list of two files"
    al_conc_groups_1 = re.split(
        r"(?:1|2)\.(fq|fastq)(?:\.gz)?$", concordant[0], flags=re.IGNORECASE
    )
    al_conc_groups_2 = re.split(
        r"(?:1|2)\.(fq|fastq)(?:\.gz)?$", concordant[1], flags=re.IGNORECASE
    )
    assert set(al_conc_groups_1) == set(
        al_conc_groups_2
    ), "input: concordant file name pair is invalid"
    al_conc_opt = (
        "al-conc-gz"
        if al_conc_groups_1[2].lower() == ".gz"
        else "al-conc-bz2" if al_conc_groups_1[2].lower() == ".bz2" else "al-conc"
    )
    al_conc_pat = f"{al_conc_groups_1[0]}%.{al_conc_groups_1[1]}"
    extra += f" --{al_conc_opt} {al_conc_pat} "

shellcmd = (
    f"(hisat2"
    f" --threads {snakemake.threads}"
    f" {reads} "
    f" -x '{index}'"
    f" {extra}"
    f" | samtools view -Sbh -o {snakemake.output[0]} -) "
    f") {log}"
)
shellcmd = re.sub(r"\s+", " ", shellcmd)
with open(snakemake.log[0], "wt") as log_fh:
    log_fh.write(f"{shellcmd}\n")

shell(shellcmd)
