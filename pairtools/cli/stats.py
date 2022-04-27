#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import sys
import click

import pandas as pd
from ..lib import fileio, pairsam_format, headerops
from . import cli, common_io_options

from ..lib.stats import PairCounter, do_merge

from .._logging import get_logger
logger = get_logger()

UTIL_NAME = "pairtools_stats"


@cli.command()
@click.argument("input_path", type=str, nargs=-1, required=False)
@click.option("-o", "--output", type=str, default="", help="output stats tsv file.")
@click.option(
    "--merge",
    is_flag=True,
    help="If specified, merge multiple input stats files instead of calculating"
    " statistics of a .pairs/.pairsam file. Merging is performed via summation of"
    " all overlapping statistics. Non-overlapping statistics are appended to"
    " the end of the file.",
)
@click.option(
    "--with-chromsizes/--no-chromsizes",
    is_flag=True,
    default=True,
    help="If specified, merge multiple input stats files instead of calculating"
    " statistics of a .pairs/.pairsam file. Merging is performed via summation of"
    " all overlapping statistics. Non-overlapping statistics are appended to"
    " the end of the file.",
)
@click.option(
    "--yaml/--no-yaml",
    is_flag=True,
    default=False,
    help="Output stats in yaml format instead of table. ",
)
@click.option(
    "--bytile-dups/--no-bytile-dups",
    default=False,
    help="If enabled, will analyse by-tile duplication statistics to estimate"
    " library complexity more accurately."
    " Requires parent_readID column to be saved by dedup (will be ignored otherwise)"
    " Saves by-tile stats into --output_bytile-stats stream, or regular output if --output_bytile-stats is not provided.",
)
@click.option(
    "--output-bytile-stats",
    default="",
    required=False,
    help="output file for tile duplicate statistics."
    " If file exists, it will be open in the append mode."
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4c-compressed."
    " By default, by-tile duplicate statistics are not printed."
    " Note that the readID and parent_readID should be provided and contain tile information for this option.",
)
@common_io_options
def stats(input_path, output, merge, bytile_dups, output_bytile_stats, **kwargs):
    """Calculate pairs statistics.

    INPUT_PATH : by default, a .pairs/.pairsam file to calculate statistics.
    If not provided, the input is read from stdin.
    If --merge is specified, then INPUT_PATH is interpreted as an arbitrary number
    of stats files to merge.

    The files with paths ending with .gz/.lz4 are decompressed by bgzip/lz4c.
    """

    stats_py(
        input_path, output, merge, bytile_dups, output_bytile_stats, **kwargs,
    )


def stats_py(
    input_path, output, merge, bytile_dups, output_bytile_stats, **kwargs
):
    if merge:
        do_merge(output, input_path, **kwargs)
        return

    if len(input_path) == 0:
        raise ValueError(f"No input paths: {input_path}")

    instream = fileio.auto_open(
        input_path[0],
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
    if bytile_dups and not output_bytile_stats:
        output_bytile_stats = outstream
    if output_bytile_stats:
        bytile_dups = True

    header, body_stream = headerops.get_header(instream)
    cols = headerops.extract_column_names(header)

    # Check necessary columns for reporting by-tile stats:
    if bytile_dups and "parent_readID" not in cols:
        logger.warning(
            "No 'parent_readID' column in the file, not generating duplicate stats."
        )
        bytile_dups = False

    # new stats class stuff would come here ...
    stats = PairCounter(bytile_dups=bytile_dups)

    # Collecting statistics
    for chunk in pd.read_table(body_stream, names=cols, chunksize=100_000):
        stats.add_pairs_from_dataframe(chunk)

    if kwargs.get("with_chromsizes", True):
        chromsizes = headerops.extract_chromsizes(header)
        stats.add_chromsizes(chromsizes)

    if bytile_dups:
        stats.save_bytile_dups(output_bytile_stats)

    # save statistics to file ...
    stats.save(outstream, yaml=kwargs.get("yaml", False))

    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()


if __name__ == "__main__":
    stats()
