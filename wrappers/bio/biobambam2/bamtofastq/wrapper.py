__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "BSD 3-Clause"

import re

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

log = snakemake.log_fmt_shell(
    stdout=True if per_readgrp else False, stderr=True, append=True
)

shellcmd = f"bamtofastq filename={snakemake.input} {output} {extra} {log}"
shellcmd = re.sub(r"\s+", " ", shellcmd)
with open(snakemake.log[0], "wt") as log_fh:
    log_fh.write(f"{shellcmd}\n")

shell(shellcmd)
