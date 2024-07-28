from contextlib import redirect_stdout, redirect_stderr
from urllib.request import urlcleanup, urlretrieve

with open(snakemake.log[0], "wt") as log_fh:
    with redirect_stdout(log_fh), redirect_stderr(log_fh):
        urlretrieve(snakemake.params[0], filename=snakemake.output[0])
        urlcleanup()
