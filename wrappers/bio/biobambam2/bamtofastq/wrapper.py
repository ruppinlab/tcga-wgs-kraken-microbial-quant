__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

per_readgrp = snakemake.params.get("per_readgrp", False)
paired_end = snakemake.params.get("paired_end", False)

extra = snakemake.params.get("extra", "")

if per_readgrp:
    outdir = snakemake.params.get("outdir")
    assert (
        outdir is not None
    ), "params: outdir is a required input parameter when per_readgrp=True"
    output = f"outputdir={outdir}"
elif paired_end:
    output = (
        f"F={snakemake.output[0]} F2={snakemake.output[1]} "
        f"O={snakemake.output[2]} 02={snakemake.output[3]}"
    )
else:
    output = f"S={snakemake.output[0]}"

log = snakemake.log_fmt_shell(stdout=True if per_readgrp else False, stderr=True)

shell("bamtofastq filename={snakemake.input} {output} {extra} {log}")
