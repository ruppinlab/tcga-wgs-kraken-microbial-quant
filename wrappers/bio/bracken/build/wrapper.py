__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "MIT"

import re

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

db = snakemake.input.get("db") or snakemake.params.get("db")
assert db is not None, "input/params: db is a required parameter"

read_length = snakemake.params.get("readlen")
assert read_length is not None, "params: readlen is a required parameter"
assert re.search(r"^\d+$", str(read_length)), "params: readlen must be an integer"

bracken_build = snakemake.params.get("bracken_build", "bracken-build")

shellcmd = (
    f"{bracken_build}"
    f" -d {db}"
    f" -k {snakemake.params.klen}"
    f" -l {read_length}"
    f" -y {snakemake.params.ktype}"
    f" -t {snakemake.threads}"
    f" {log}"
)
shellcmd = re.sub(r"\s+", " ", shellcmd)
with open(snakemake.log[0], "wt") as log_fh:
    log_fh.write(f"{shellcmd}\n")
shell(shellcmd)
