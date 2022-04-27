#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import sys
import click
import warnings

import pandas as pd
from ..lib import fileio, pairsam_format, headerops
from . import cli, common_io_options

from ..lib.stats import PairCounter, do_merge

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
    "--bytile-dups/--no-bytile-dups",
    default=False,
    help="If enabled, will analyse by-tile duplication statistics to estimate"
    " library complexity more accurately."
    " Requires parent_readID column to be saved by dedup (will be ignored otherwise)",
)
@click.option(
    "--output-bytile-dups-stats",
    type=str,
    default="",
    help="If specified, will analyse by-tile duplication statistics to estimate"
    " library complexity more accurately and will save details to this path."
    " Requires parent_readID column to be saved by dedup (will be ignored otherwise)",
)
@common_io_options
def stats(input_path, output, merge, bytile_dups, output_bytile_dups_stats, **kwargs):
    """Calculate pairs statistics.

    INPUT_PATH : by default, a .pairs/.pairsam file to calculate statistics.
    If not provided, the input is read from stdin.
    If --merge is specified, then INPUT_PATH is interpreted as an arbitrary number
    of stats files to merge.

    The files with paths ending with .gz/.lz4 are decompressed by bgzip/lz4c.
    """
    stats_py(
        input_path, output, merge, bytile_dups, output_bytile_dups_stats, **kwargs,
    )


def stats_py(
    input_path, output, merge, bytile_dups, output_bytile_dups_stats, **kwargs
):
    if merge:
        do_merge(output, input_path, **kwargs)
        return

    instream = (
        fileio.auto_open(
            input_path[0],
            mode="r",
            nproc=kwargs.get("nproc_in"),
            command=kwargs.get("cmd_in", None),
        )
        if input_path
        else sys.stdin
    )
    outstream = (
        fileio.auto_open(
            output,
            mode="w",
            nproc=kwargs.get("nproc_out"),
            command=kwargs.get("cmd_out", None),
        )
        if output
        else sys.stdout
    )

    header, body_stream = headerops.get_header(instream)
    cols = headerops.extract_column_names(header)

    if (bytile_dups or output_bytile_dups_stats) and "parent_readID" not in cols:
        warnings.warn(
            "No 'parent_readID' column in the file, not generating duplicate stats."
        )
        bytile_dups = False
        output_bytile_dups_stats = False
    # new stats class stuff would come here ...
    stats = PairCounter()

    # Collecting statistics
    for chunk in pd.read_table(body_stream, names=cols, chunksize=100_000):
        stats.add_pairs_from_dataframe(chunk)

    if output_bytile_dups_stats:
        stats.save_bytile_dups()
    stats.calculate_summaries()

    # save statistics to file ...
    stats.save(outstream)

    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()


if __name__ == "__main__":
    stats()
