__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "MIT"

import re
from os.path import exists, join
from tempfile import gettempdir, NamedTemporaryFile, TemporaryDirectory
from shutil import move

from snakemake.shell import shell
from snakemake.utils import makedirs

rg_meta_df = snakemake.params.get("rg_meta_df")
assert rg_meta_df is not None, "params: rg_meta_df is a required parameter"

per_readgrp = (
    rg_meta_df["is_paired_end"].nunique() > 1 or rg_meta_df["read_length"].nunique() > 1
)
paired_end = rg_meta_df["is_paired_end"].all()

extra = snakemake.params.get("extra", "")

if per_readgrp:
    suffixes = snakemake.params.get("suffixes")
    assert (
        suffixes is not None
    ), "params: suffixes is a required parameter when per readgroup output"
    makedirs(snakemake.output[0])
    extra = f"collate=1 combs=1 {extra} outputperreadgroup=1"
    output = f"outputdir={snakemake.output[0]}"
    for suffix_opt in suffixes:
        output += f" outputperreadgroupsuffix{suffix_opt}={suffixes[suffix_opt]}"
elif paired_end:
    extra = f"collate=1 combs=1 {extra}"
    output = (
        f"F={snakemake.output.F} F2={snakemake.output.F2} "
        f"O={snakemake.params.O} 02={snakemake.params.O2} "
        f"> /dev/null"
    )
else:
    output = f"> {snakemake.output[0]}"

log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)

with TemporaryDirectory(dir=snakemake.resources.get("tmpdir", gettempdir())) as tmpdir:
    with NamedTemporaryFile(
        dir=tmpdir, prefix="bamtofastq_", delete=False, delete_on_close=False
    ) as tmpfile:
        shellcmd = (
            f"bamtofastq"
            f" filename={snakemake.input[0]}"
            f" {extra}"
            f" T={tmpfile.name}"
            f" {output}"
            f" {log}"
        )
        shellcmd = re.sub(r"\s+", " ", shellcmd)
        with open(snakemake.log[0], "wt") as log_fh:
            log_fh.write(f"{shellcmd}\n")
        shell(shellcmd)

if per_readgrp:
    with open(snakemake.log[0], "at") as log_fh:
        for _, rg in rg_meta_df.iterrows():
            for suffix_opt in suffixes:
                outfile_src = join(
                    snakemake.output[0],
                    f"{rg['read_group_name']}{suffixes[suffix_opt]}",
                )
                outfile_dst = join(
                    snakemake.output[0], f"{rg['read_group_id']}{suffixes[suffix_opt]}"
                )
                if exists(outfile_src):
                    move(outfile_src, outfile_dst)
                    log_fh.write(f"{outfile_src} -> {outfile_dst}\n")
