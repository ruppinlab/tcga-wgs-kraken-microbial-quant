__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "MIT"

import re

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

db = snakemake.input.get("db") or snakemake.params.get("db")
assert db is not None, "input/params: db is a required parameter"

read_length_file = snakemake.input.get("readlen")
if read_length_file is not None:
    with open(read_length_file, "r") as fh:
        read_length = re.sub(r"\D+", "", fh.readline())
else:
    read_length = snakemake.params.get("readlen")

assert read_length is not None, "input/params: readlen is a required parameter"
assert re.search(r"^\d+$", str(read_length)), "params: readlen must be an integer"

db_read_lengths = snakemake.params.get("db_readlens")
assert db_read_lengths is not None, "params: db_readlens is a required parameter"

# find closest db read length <= sample read length
read_length = max([l for l in db_read_lengths if int(read_length) >= int(l)])

level = snakemake.params.get("level", "S")
read_thres = snakemake.params.get("read_thres", 10)
kmer_thres = snakemake.params.get("kmer_thres", 0)

bracken = snakemake.params.get("bracken", "bracken")

shellcmd = (
    f"{bracken}"
    f" -d {db}"
    f" -i {snakemake.input.report}"
    f" -o {snakemake.output.counts}"
    f" -w {snakemake.output.report}"
    f" -r {read_length}"
    f" -l {level}"
    f" -t {read_thres}"
    f" -m {kmer_thres}"
    f" {log}"
)
shellcmd = re.sub(r"\s+", " ", shellcmd)
with open(snakemake.log[0], "wt") as log_fh:
    log_fh.write(f"{shellcmd}\n")
shell(shellcmd)
