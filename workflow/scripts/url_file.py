from os import system
from shutil import copyfileobj
from urllib.error import URLError
from urllib.request import Request, urlopen

req = Request(
    snakemake.params.url,
    data=None,
    headers=snakemake.params.headers,
    method=snakemake.params.method,
)

with urlopen(req) as res, open(snakemake.output[0], "wb") as out_file:
    copyfileobj(res, out_file)

if system(f"samtools quickcheck -v '{snakemake.output[0]}'") != 0:
    raise URLError("BAM download incomplete, quickcheck failed")
