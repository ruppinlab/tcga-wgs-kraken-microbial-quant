__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "MIT"

import re

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

extra = snakemake.params.get("extra", "")

combine_kreports = snakemake.params.get("combine_kreports", "combine_kreports.py")

shellcmd = (
    f"{combine_kreports}"
    f" {extra}"
    f" --report-files {snakemake.input}"
    f" --output {snakemake.output}"
    f" {log}"
)
shellcmd = re.sub(r"\s+", " ", shellcmd)
with open(snakemake.log[0], "wt") as log_fh:
    log_fh.write(f"{shellcmd}\n")
shell(shellcmd)
