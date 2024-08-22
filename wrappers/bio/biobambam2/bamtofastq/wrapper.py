__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "MIT"

import re
from os import remove
from os.path import exists, join
from tempfile import gettempdir, NamedTemporaryFile, TemporaryDirectory
from shutil import move

from snakemake.shell import shell
from snakemake.utils import makedirs

per_readgrp = snakemake.params.get("per_readgrp", False)
paired_end = snakemake.params.get("paired_end", False)
extra = snakemake.params.get("extra", "")

if per_readgrp:
    readgrp_names = snakemake.params.get("readgrp_names")
    assert (
        readgrp_names is not None
    ), "params: readgrp_names is a required parameter when per_readgrp=True"
    readgrp_ids = snakemake.params.get("readgrp_ids")
    assert (
        readgrp_ids is not None
    ), "params: readgrp_ids is a required parameter when per_readgrp=True"
    assert (
        isinstance(readgrp_names, (list, tuple))
        and isinstance(readgrp_ids, (list, tuple))
        and len(readgrp_names) == len(readgrp_ids)
    ), "params: readgrp_names and readgrp_ids must be lists of same length"
    excl_readgrp_ids = snakemake.params.get("excl_readgrp_ids", [])
    assert isinstance(
        excl_readgrp_ids, (list, tuple)
    ), "params: excl_readgrp_ids must be list-like"
    suffixes = snakemake.params.get("suffixes")
    assert (
        suffixes is not None
    ), "params: suffixes is a required parameter when per_readgrp=True"
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
        for rg_name, rg_id in zip(readgrp_names, readgrp_ids):
            for suffix_opt in suffixes:
                outfile_src = join(
                    snakemake.output[0], f"{rg_name}{suffixes[suffix_opt]}"
                )
                outfile_dst = join(
                    snakemake.output[0], f"{rg_id}{suffixes[suffix_opt]}"
                )
                if exists(outfile_src):
                    move(outfile_src, outfile_dst)
                    log_fh.write(f"{outfile_src} -> {outfile_dst}\n")
                    if excl_readgrp_ids and rg_id in excl_readgrp_ids:
                        remove(outfile_dst)
                        log_fh.write(f"Removed/excluded {outfile_dst}\n")
