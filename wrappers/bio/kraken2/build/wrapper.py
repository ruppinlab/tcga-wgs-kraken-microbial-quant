__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

taskopt = snakemake.params.get("taskopt")
assert taskopt is not None, "params: taskopt is a required parameter"
assert taskopt in (
    "--download-taxonomy",
    "--download-library",
), "params: invalid taskopt"
if taskopt == "--download-library":
    lib = snakemake.params.get("lib")
    assert (
        lib is not None
    ), "params: lib is a required parameter when taskopt=--download-library"
    taskopt = f"{taskopt} {lib}"

extra = snakemake.params.get("extra", "")

shell(
    "kraken2-build"
    " {taskopt}"
    " --db {snakemake.input.db}"
    " --threads {snakemake.threads}"
    " {extra}"
    " {log}"
)
