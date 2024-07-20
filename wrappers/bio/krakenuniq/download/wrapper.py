__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

db = snakemake.params.get("db")
assert db is not None, "params: db is a required parameter"

lib = snakemake.params.get("lib")
assert lib is not None, "params: lib is a required parameter"

extra = snakemake.params.get("extra", "")

shell(
    "krakenuniq-download"
    " --db {db}"
    " --threads {snakemake.threads}"
    " {extra}"
    " {lib}"
    " {log}"
)
