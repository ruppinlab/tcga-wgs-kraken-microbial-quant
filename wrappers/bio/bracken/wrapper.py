__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")

shell(
    "bracken"
    " -d {snakemake.input.db}"
    " -i {snakemake.input.report}"
    " -o {snakemake.output}"
    " -r {snakemake.params.readlen}"
    " -l {snakemake.params.level}"
    " -t {snakemake.params.threshold}"
    " {extra}"
    " {log}"
)
