__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "MIT"

import re

from snakemake.shell import shell

fastq_ext_regex = re.compile(r"_(?:1|2)(\.(?:fq|fastq))(?:\.gz)?$", flags=re.IGNORECASE)

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

paired_end = snakemake.params.get("paired_end", False)

output = snakemake.output.get("output")
if output is None:
    output = snakemake.params.get("output")
assert output is not None, "output/params: output is a required named parameter"
output = f"--output {output}"
classif = snakemake.output.get("classif")
if classif is not None:
    if paired_end:
        classif = re.sub(fastq_ext_regex, r"#\1", classif[0], count=1)
    output += f" --classified-out {classif}"
unclassif = snakemake.output.get("unclassif")
if unclassif is not None:
    if paired_end:
        unclassif = re.sub(fastq_ext_regex, r"#\1", unclassif[0], count=1)
    output += f" --unclassified-out {unclassif}"
report = snakemake.output.get("report")
if report is not None:
    report = f"--report {report}"

extra = snakemake.params.get("extra", "")
if paired_end:
    extra = f"--paired {extra}"

shellcmd = (
    f"kraken2"
    f" --db {snakemake.params.db}"
    f" --threads {snakemake.threads}"
    f" {output}"
    f" {report}"
    f" {extra}"
    f" {snakemake.input.fqs}"
    f" {log}"
)
shellcmd = re.sub(r"\s+", " ", shellcmd)
with open(snakemake.log[0], "wt") as log_fh:
    log_fh.write(f"{shellcmd}\n")

shell(shellcmd)
