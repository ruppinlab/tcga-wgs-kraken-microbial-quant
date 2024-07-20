__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "BSD 3-Clause"

import re

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

read_length = snakemake.params.get("readlen")
if read_length is None:
    read_length_file = snakemake.input.get("readlen")
    assert read_length_file is not None, "input/params: readlen is a required parameter"
    with open(read_length_file, "r") as fh:
        read_length = re.sub("\D+", "", fh.readline())

db_read_lengths = snakemake.params.get("db_readlens")
assert db_read_lengths is not None, "params: db_readlens is a required parameter"

# find closest db read length to sample read length
read_length = min(db_read_lengths, lambda x: abs(x - read_length))

level = snakemake.params.get("level", "S")
threshold = snakemake.params.get("threshold", 0.0)

shell(
    "bracken"
    " -d {snakemake.input.db}"
    " -i {snakemake.input.report}"
    " -o {snakemake.output.counts}"
    " -w {snakemake.output.report}"
    " -r {read_length}"
    " -l {level}"
    " -t {threshold}"
    " {log}"
)
