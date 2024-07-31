from shutil import copyfileobj
from urllib.request import Request, urlopen

req = Request(
    f"https://api.gdc.cancer.gov/slicing/view/{snakemake.params['bam_id']}?region=unmapped",
    data=None,
    headers={"X-Auth-Token": snakemake.params["token"]},
    method="GET",
)

with urlopen(req) as res, open(snakemake.output[0], "wb") as out_file:
    copyfileobj(res, out_file)
