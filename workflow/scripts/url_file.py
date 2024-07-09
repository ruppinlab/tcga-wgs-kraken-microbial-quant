from urllib.request import urlcleanup, urlretrieve

urlretrieve(snakemake.params[0], filename=snakemake.output[0])
urlcleanup()
