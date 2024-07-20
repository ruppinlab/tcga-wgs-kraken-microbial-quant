__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "BSD 3-Clause"

import re

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

db = snakemake.input.get("db") or snakemake.input[0]
assert db is not None, "input: db is a required position 0 or named parameter"

read_length = snakemake.params.get("readlen")
if read_length is None:
    read_length_file = snakemake.input.get("readlen")
    assert read_length_file is not None, "input/params: readlen is a required parameter"
    with open(read_length_file, "r") as fh:
        read_length = re.sub("\D+", "", fh.readline())

shell(
    "bracken-build"
    " -d {db}"
    " -k {snakemake.params.klen}"
    " -l {read_length}"
    " -y {snakemake.params.ktype}"
    " -t {snakemake.threads}"
    " {log}"
)
