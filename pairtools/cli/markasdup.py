#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import click

from ..lib import fileio, pairsam_format, headerops
from . import cli, common_io_options

from ..lib.dedup import mark_split_pair_as_dup


UTIL_NAME = "pairtools_markasdup"


@cli.command()
@click.argument("pairsam_path", type=str, required=False)
@click.option(
    "-o",
    "--output",
    type=str,
    default="",
    help="output .pairsam file."
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4c-compressed."
    " By default, the output is printed into stdout.",
)
@common_io_options
def markasdup(pairsam_path, output, **kwargs):
    """Tag all pairs in the input file as duplicates.

    Change the type of all pairs inside a .pairs/.pairsam file to DD. If sam
    entries are present, change the pair type in the Yt SAM tag to 'Yt:Z:DD'.

    PAIRSAM_PATH : input .pairs/.pairsam file. If the path ends with .gz, the
    input is gzip-decompressed. By default, the input is read from stdin.
    """
    markasdup_py(pairsam_path, output, **kwargs)


def markasdup_py(pairsam_path, output, **kwargs):
    instream = fileio.auto_open(
        pairsam_path,
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
    outstream.writelines((l + "\n" for l in header))

    for line in body_stream:
        cols = line.rstrip('\n').split(pairsam_format.PAIRSAM_SEP)
        mark_split_pair_as_dup(cols)

        outstream.write(pairsam_format.PAIRSAM_SEP.join(cols))
        outstream.write("\n")

    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()


if __name__ == "__main__":
    markasdup()
