__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

outdir = (
    snakemake.output.get("outdir")
    if snakemake.output.get("outdir")
    else snakemake.params.get("outdir")
)

extra = snakemake.params.get("extra", "")

shell(
    "bamtofastq"
    " filename={snakemake.input}"
    " outputdir={outdir}"
    " {extra}"
    " {log}"
)
