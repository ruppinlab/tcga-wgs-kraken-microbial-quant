__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")
report = (
    f"--report-file {snakemake.output.get("report")}"
    if snakemake.output.get("report")
    else ""
)

shell(
    "kraken2"
    " --db {snakemake.input.db}"
    " --threads {snakemake.threads}"
    " --output {output}"
    " {extra}"
    " {report}"
    " {snakemake.input.fqs}"
    " {log}"
)
