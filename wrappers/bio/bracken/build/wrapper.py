__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "MIT"

import re

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

db = snakemake.input.get("db") or snakemake.params.get("db")
assert db is not None, "input/params: db is a required parameter"

ktype = snakemake.params.get("ktype", "kraken2")
assert ktype in (
    "kraken",
    "kraken2",
    "krakenuniq",
), "params: ktype invalid must = kraken | kraken2 | krakenuniq"

klen = snakemake.params.get("klen", 35 if ktype == "kraken2" else 31)

read_length = snakemake.params.get("readlen", 100)
assert re.search(r"^\d+$", str(read_length)), "params: readlen must be an integer"

db_only = snakemake.params.get("db_only", False)
db_only = "-o" if db_only else ""

bracken_build = snakemake.params.get("bracken_build", "bracken-build")

shellcmd = (
    f"{bracken_build}"
    f" -d {db}"
    f" -k {klen}"
    f" -l {read_length}"
    f" -y {ktype}"
    f" -t {snakemake.threads}"
    f" {db_only}"
    f" {log}"
)
shellcmd = re.sub(r"\s+", " ", shellcmd)
with open(snakemake.log[0], "wt") as log_fh:
    log_fh.write(f"{shellcmd}\n")
shell(shellcmd)
