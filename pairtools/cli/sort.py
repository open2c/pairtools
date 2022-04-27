#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import sys
import click
import subprocess
import shutil
import warnings

from ..lib import fileio, pairsam_format, headerops
from . import cli, common_io_options

UTIL_NAME = "pairtools_sort"


@cli.command()
@click.argument("pairs_path", type=str, required=False)
@click.option(
    "-o",
    "--output",
    type=str,
    default="",
    help="output pairs file."
    " If the path ends with .gz or .lz4, the output is compressed by bgzip "
    "or lz4, correspondingly. By default, the output is printed into stdout.",
)
@click.option(
    "--nproc",
    type=int,
    default=8,
    show_default=True,
    help="Number of processes to split the sorting work between.",
)
@click.option(
    "--tmpdir",
    type=str,
    default="",
    help="Custom temporary folder for sorting intermediates.",
)
@click.option(
    "--memory",
    type=str,
    default="2G",
    show_default=True,
    help="The amount of memory used by default.",
)
@click.option(
    "--compress-program",
    type=str,
    default="auto",
    show_default=True,
    help="A binary to compress temporary sorted chunks. "
    "Must decompress input when the flag -d is provided. "
    "Suggested alternatives: gzip, lzop, lz4c, snzip. "
    'If "auto", then use lz4c if available, and gzip '
    "otherwise.",
)
@common_io_options
def sort(pairs_path, output, nproc, tmpdir, memory, compress_program, **kwargs):
    """Sort a .pairs/.pairsam file.

    Sort pairs in the lexicographic order along chrom1 and chrom2, in the
    numeric order along pos1 and pos2 and in the lexicographic order along
    pair_type.

    PAIRS_PATH : input .pairs/.pairsam file. If the path ends with .gz or .lz4, the
    input is decompressed by bgzip or lz4c, correspondingly. By default, the
    input is read as text from stdin.
    """
    sort_py(pairs_path, output, nproc, tmpdir, memory, compress_program, **kwargs)


def sort_py(pairs_path, output, nproc, tmpdir, memory, compress_program, **kwargs):

    instream = fileio.auto_open(
        pairs_path,
        mode="r",
        nproc=kwargs.get("nproc_in"),
        command=kwargs.get("cmd_in", None),
    )
    outstream = fileio.auto_open(
        output,
        mode="w",
        nproc=kwargs.get("nproc_out"),
        command=kwargs.get("cmd_out", None),
    )

    header, body_stream = headerops.get_header(instream)
    header = headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    header = headerops.mark_header_as_sorted(header)

    outstream.writelines((l + "\n" for l in header))

    outstream.flush()

    if compress_program == "auto":
        if shutil.which("lz4c") is not None:
            compress_program = "lz4c"
        else:
            warnings.warn(
                "lz4c is not found. Using gzip for compression of sorted chunks, "
                "which results in a minor decrease in performance. Please install "
                "lz4c for faster sorting."
            )
            compress_program = "gzip"

    command = r"""
        /bin/bash -c 'export LC_COLLATE=C; export LANG=C; sort 
        -k {0},{0} -k {1},{1} -k {2},{2}n -k {3},{3}n -k {4},{4} 
        --stable
        --field-separator=$'\''{5}'\'' 
        {6}
        {7}
        -S {8}
        {9}
        """.replace(
        "\n", " "
    ).format(
        pairsam_format.COL_C1 + 1,
        pairsam_format.COL_C2 + 1,
        pairsam_format.COL_P1 + 1,
        pairsam_format.COL_P2 + 1,
        pairsam_format.COL_PTYPE + 1,
        pairsam_format.PAIRSAM_SEP_ESCAPE,
        " --parallel={} ".format(nproc) if nproc > 0 else " ",
        " --temporary-directory={} ".format(tmpdir) if tmpdir else " ",
        memory,
        (
            " --compress-program={} ".format(compress_program)
            if compress_program
            else " "
        ),
    )
    command += "'"

    with subprocess.Popen(
        command, stdin=subprocess.PIPE, bufsize=-1, shell=True, stdout=outstream
    ) as process:
        stdin_wrapper = io.TextIOWrapper(process.stdin, "utf-8")
        for line in body_stream:
            stdin_wrapper.write(line)
        stdin_wrapper.flush()
        process.communicate()

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()


if __name__ == "__main__":
    sort()
