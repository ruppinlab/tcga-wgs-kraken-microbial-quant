__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "BSD 3-Clause"

import re

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

db = snakemake.params.get("db")
assert db is not None, "params: db is a required parameter"

lib = snakemake.params.get("lib")
assert lib is not None, "params: lib is a required parameter"

extra = snakemake.params.get("extra", "")
rsync = snakemake.params.get("rsync")
if rsync:
    extra = f"--rsync {extra}"

shellcmd = (
    f"krakenuniq-download"
    f" --db {db}"
    f" --threads {snakemake.threads}"
    f" {extra}"
    f" {lib}"
    f" {log}"
)
shellcmd = re.sub(r"\s+", " ", shellcmd)
with open(snakemake.log[0], "wt") as log_fh:
    log_fh.write(f"{shellcmd}\n")

shell(shellcmd)
