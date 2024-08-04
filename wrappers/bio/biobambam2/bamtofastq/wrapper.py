__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "MIT"

import re
from os.path import exists, isdir, join
from tempfile import gettempdir, NamedTemporaryFile, TemporaryDirectory
from shutil import move

from snakemake.shell import shell

per_readgrp = snakemake.params.get("per_readgrp", False)
paired_end = snakemake.params.get("paired_end", False)

extra = snakemake.params.get("extra", "")

if per_readgrp:
    outputdir = snakemake.params.get("outputdir")
    assert (
        outputdir is not None
    ), "params: outputdir is a required parameter when per_readgrp=True"
    assert isdir(outputdir), "params: outputdir must be an existing directory"
    output = f"outputdir={outputdir}"
    rg_meta_df = snakemake.params.get("rg_meta_df")
    assert (
        rg_meta_df is not None
    ), "params: rg_meta_df is a required parameter when per_readgrp=True"
    suffixes = {}
    suffix_opts = (
        [
            "outputperreadgroupsuffixF",
            "outputperreadgroupsuffixF2",
            "outputperreadgroupsuffixO",
            "outputperreadgroupsuffixO2",
        ]
        if paired_end
        else ["outputperreadgroupsuffixS"]
    )
    for suffix_opt in suffix_opts:
        suffix = snakemake.params.get(suffix_opt)
        assert (
            suffix is not None
        ), f"params: {suffix_opt} is a required parameter when per_readgrp=True"
        suffixes[suffix_opt] = suffix
        extra = f" {suffix_opt}={suffix} {extra}"
    if paired_end:
        extra = f" collate=1 {extra}"
elif paired_end:
    output = (
        f"F={snakemake.output.F} F2={snakemake.output.F2} "
        f"O={snakemake.params.O} 02={snakemake.params.O2} "
        f"> /dev/null"
    )
    extra = f" collate=1 {extra}"
else:
    output = f"> {snakemake.output[0]}"

log = snakemake.log_fmt_shell(
    stdout=False if paired_end or per_readgrp else True, stderr=True, append=True
)

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

with open(snakemake.log[0], "at") as log_fh:
    if per_readgrp:
        for _, rg in rg_meta_df.iterrows():
            for suffix_opt in suffixes.keys():
                outfile_src = join(
                    outputdir, f"{rg['read_group_name']}{suffixes[suffix_opt]}"
                )
                outfile_dst = join(
                    outputdir, f"{rg['read_group_id']}{suffixes[suffix_opt]}"
                )
                if exists(outfile_src):
                    move(outfile_src, outfile_dst)
                    log_fh.write(f"Moved {outfile_src} -> {outfile_dst}\n")
