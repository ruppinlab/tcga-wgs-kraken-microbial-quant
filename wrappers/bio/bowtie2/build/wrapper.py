__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "BSD 3-Clause"

import re
from os.path import commonprefix

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

extra = snakemake.params.get("extra", "")

index = commonprefix(snakemake.output).rstrip(".")

shellcmd = (
    f"bowtie2-build"
    f" --threads {snakemake.threads}"
    f" {extra}"
    f" {snakemake.input.ref}"
    f" {index}"
    f" {log}"
)
shellcmd = re.sub(r"\s+", " ", shellcmd)
with open(snakemake.log[0], "wt") as log_fh:
    log_fh.write(f"{shellcmd}\n")

shell(shellcmd)
