__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "MIT"

import re

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

paired_end = snakemake.params.get("paired_end", False)

output = snakemake.output.get("output")
if output is None:
    output = snakemake.params.get("output")
assert output is not None, "output/params: output is a required named parameter"
output = f"--output {output}"
classif = snakemake.output.get("classif")
if classif is not None:
    output += f" --classified-out {classif}"
unclassif = snakemake.output.get("unclassif")
if unclassif is not None:
    output += f" --unclassified-out {unclassif}"
report = snakemake.output.get("report")
if report is not None:
    report = f"--report-file {report}"

extra = snakemake.params.get("extra", "")
if paired_end:
    extra = f"--paired {extra}"

shellcmd = (
    f"krakenuniq"
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
