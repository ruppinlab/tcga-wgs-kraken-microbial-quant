__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

db = snakemake.params.get("db")
assert db is not None, "params: db is a required parameter"

taskopt = snakemake.params.get("taskopt")
assert taskopt is not None, "params: taskopt is a required parameter"
assert taskopt in (
    "--download-taxonomy",
    "--download-library",
    "--build",
), "params: invalid taskopt"
if taskopt == "--download-library":
    lib = snakemake.params.get("lib")
    assert (
        lib is not None
    ), "params: lib is a required parameter when taskopt=--download-library"
    taskopt = f"{taskopt} {lib}"

extra = snakemake.params.get("extra", "")
protein = snakemake.params.get("protein")
if protein is not None:
    extra = f"--protein {extra}"

# workaround for Kraken2 issue with --download-library human
kraken2_build = (
    "yes y | kraken2-build"
    if taskopt == "--download-library" and lib == "human"
    else "kraken2-build"
)

shell(
    "{kraken2_build}"
    " {taskopt}"
    " --db {db}"
    " --threads {snakemake.threads}"
    " {extra}"
    " {log}"
)
