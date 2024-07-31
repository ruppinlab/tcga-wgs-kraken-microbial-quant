__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "MIT"

import re

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

fq1 = snakemake.input.get("fq1") or snakemake.input[0]
assert fq1 is not None, "input: fq1 is a required named or positional input parameter"

fq2 = snakemake.input.get("fq2") or snakemake.input[1]
if fq2:
    in2 = f"in2='{fq2}'"
else:
    in2 = ""

shellcmd = f"readlength.sh in='{fq1}' {in2} bin=1 out='{snakemake.output}' {log}"
shellcmd = re.sub(r"\s+", " ", shellcmd)
with open(snakemake.log[0], "wt") as log_fh:
    log_fh.write(f"{shellcmd}\n")

shell(shellcmd)
