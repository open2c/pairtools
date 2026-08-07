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
    " the end of the file. Supported for tsv stats with single filter.",
)
@click.option(
    "--single-mapped-by-side",
    is_flag=True,
    help="If specified, count single-mapped reads separately for read1 and read2.",
)
@click.option(
    "--n-dist-bins-decade",
    type=int,
    default=PairCounter.N_DIST_BINS_DECADE_DEFAULT,
    show_default=True,
    required=False,
    help="Number of bins to split the distance range in log10-space, specified per a factor of 10 difference.",
)
@click.option(
    "--with-chromsizes/--no-chromsizes",
    is_flag=True,
    default=True,
    help="If enabled, will store sizes of chromosomes from the header of the pairs file"
    " in the stats file.",
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
# Filtering options:
@click.option(
    "--filter",
    default=None,
    required=False,
    multiple=True,
    help="Filters with conditions to apply to the data (similar to `pairtools select`). "
    "For non-YAML output only the first filter will be reported. "
    """Example: pairtools stats --yaml --filter 'unique:(pair_type=="UU")' --filter 'close:(pair_type=="UU") and (abs(pos1-pos2)<10)' test.pairs """,
)
@click.option(
    "--engine",
    default="pandas",
    required=False,
    help="Engine for regular expression parsing. "
    "Python will provide you regex functionality, while pandas does not accept custom funtctions and works faster. ",
)
@click.option(
    "--chrom-subset",
    type=str,
    default=None,
    required=False,
    help="A path to a chromosomes file (tab-separated, 1st column contains "
    "chromosome names) containing a chromosome subset of interest. "
    "If provided, additionally filter pairs with both sides originating from "
    "the provided subset of chromosomes. This operation modifies the #chromosomes: "
    "and #chromsize: header fields accordingly.",
)
@click.option(
    "--startup-code",
    type=str,
    default=None,
    required=False,
    help="An auxiliary code to execute before filtering. "
    "Use to define functions that can be evaluated in the CONDITION statement",
)
@click.option(
    "-t",
    "--type-cast",
    type=(str, str),
    default=(),
    multiple=True,
    help="Cast a given column to a given type. By default, only pos and mapq "
    "are cast to int, other columns are kept as str. Provide as "
    "-t <column_name> <type>, e.g. -t read_len1 int. Multiple entries are allowed.",
)
@common_io_options
def stats(
    input_path,
    output,
    merge,
    single_mapped_by_side,
    n_dist_bins_decade,
    bytile_dups,
    output_bytile_stats,
    filter,
    **kwargs,
):
    """Calculate pairs statistics.

    INPUT_PATH : by default, a .pairs/.pairsam file to calculate statistics.
    If not provided, the input is read from stdin.
    If --merge is specified, then INPUT_PATH is interpreted as an arbitrary number
    of stats files to merge.

    The files with paths ending with .gz/.lz4 are decompressed by bgzip/lz4c.
    """

    stats_py(
        input_path,
        output,
        merge,
        single_mapped_by_side,
        n_dist_bins_decade,
        bytile_dups,
        output_bytile_stats,
        filter,
        **kwargs,
    )


def stats_py(
    input_path,
    output,
    merge,
    single_mapped_by_side,
    n_dist_bins_decade,
    bytile_dups,
    output_bytile_stats,
    filter,
    **kwargs,
):
    if merge:
        do_merge(output, input_path, n_dist_bins_decade=n_dist_bins_decade, **kwargs)
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

    # Define filters and their properties
    first_filter_name = "no_filter"  # default filter name for full output
    if filter is not None and len(filter) > 0:
        first_filter_name = filter[0].split(":", 1)[0]
        if len(filter) > 1 and not kwargs.get("yaml", False):
            logger.warn(
                f"Output the first filter only in non-YAML output: {first_filter_name}"
            )

        filter = dict([f.split(":", 1) for f in filter])
    else:
        filter = None

    stats = PairCounter(
        single_mapped_by_side=single_mapped_by_side,
        n_dist_bins_decade=n_dist_bins_decade,
        bytile_dups=bytile_dups,
        filters=filter,
        startup_code=kwargs.get("startup_code", ""),  # for evaluation of filters
        type_cast=kwargs.get("type_cast", ()),  # for evaluation of filters
        engine=kwargs.get("engine", "pandas"),
    )

    # Collecting statistics
    for chunk in pd.read_table(body_stream, names=cols, chunksize=100_000):
        stats.add_pairs_from_dataframe(chunk)

    if kwargs.get("with_chromsizes", True):
        chromsizes = headerops.extract_chromsizes(header)
        stats.add_chromsizes(chromsizes)

    if bytile_dups:
        stats.save_bytile_dups(output_bytile_stats)

    # save statistics to file ...
    stats.save(
        outstream,
        yaml=kwargs.get("yaml", False),  # format as yaml
        filter=(
            first_filter_name if not kwargs.get("yaml", False) else None
        ),  # output only the first filter if non-YAML output
    )

    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()


if __name__ == "__main__":
    stats()
