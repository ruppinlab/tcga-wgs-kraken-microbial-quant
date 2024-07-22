__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

output = snakemake.output.get("output")
if output is None:
    output = snakemake.params.get("output")
assert output is not None, "output/params: output is a required named parameter"
output = f"--output {output}"
classif = snakemake.output.get("classif")
if classif is not None:
    output += f"--classified-out {classif}"
unclassif = snakemake.output.get("unclassif")
if unclassif is not None:
    output += f"--unclassified-out {unclassif}"
report = snakemake.output.get("report")
if report is not None:
    report = f"--report-file {report}"

extra = snakemake.params.get("extra", "")
paired_end = snakemake.params.get("paired_end", False)
if paired_end:
    extra = f"--paired {extra}"

shell(
    "kraken2"
    " --db {snakemake.input.db}"
    " --threads {snakemake.threads}"
    " {output}"
    " {report}"
    " {extra}"
    " {snakemake.input.fqs}"
    " {log}"
)
