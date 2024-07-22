__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")

shell(
    "combine_kreports.py"
    " --report-files {snakemake.input}"
    " --output {snakemake.output}"
    " {extra}"
    " {log}"
)
