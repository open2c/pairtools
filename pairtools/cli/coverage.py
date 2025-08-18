#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Command-line interface for coverage analysis.

This module provides a Click-based CLI that wraps the coverage_api module.
"""

import click
import bioframe

from ..lib.coverage import calculate_coverage, read_pairs
from . import cli, common_io_options


UTIL_NAME = "pairtools_coverage"


@cli.command()
@click.argument(
    "pairs_path",
    type=str,
    required=False,
)
@click.option(
    "--side",
    "-s",
    type=click.Choice(["both", "1", "2"]),
    default="1",
    show_default=True,
    help="both: both sides, 1: first side, 2: second side",
)
@click.option(
    "--end",
    type=click.Choice(["5", "3", "0"]),
    default="0",
    show_default=True,
    help="5', 3' end, or 0 (default) - whatever is reported as pos1 and pos2 in pairs",
)
@click.option(
    "--shift",
    type=int,
    default=0,
    show_default=True,
    help="Shift value for strand-specific adjustment. Positive values shift downstream "
    "relative to the corresponding strand, negative values shift upstream. "
    "Can be useful to e.g. account for nucleosome size in micro-C by shifting by 73 bp",
)
@click.option(
    "--window-size",
    type=int,
    default=10,
    show_default=True,
    help="Window size for binning",
)
@click.option(
    "--chunksize",
    type=int,
    default=1000000,
    show_default=True,
    help="Number of pairs to process at once",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=1,
    show_default=True,
    help="Number of threads",
)
@click.option(
    "--chromsizes",
    type=str,
    default=None,
    help="Chromosome order used to filter pairs: "
    "path to a chromosomes file (e.g. UCSC chrom.sizes or similar) whose "
    "first column lists scaffold names. If not provided, will be extracted "
    "from the header of the input pairs file.",
)
@click.option(
    "-o",
    "--output",
    type=str,
    default="",
    help="output file."
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4c-compressed."
    " By default, the output is printed into stdout.",
)
@click.option(
    "--output-bigwig",
    type=str,
    default=None,
    help="Output BigWig file (optional)",
)
@common_io_options
def coverage(
    pairs_path,
    side,
    end,
    shift,
    window_size,
    chunksize,
    threads,
    chromsizes,
    output,
    output_bigwig,
    **kwargs,
):
    """Generate coverage from pairs file.

    PAIRS_PATH : input .pairs/.pairsam file. If the path ends with .gz or .lz4, the
    input is decompressed by bgzip/lz4c.
    By default, the input is read from stdin.
    """
    coverage_py(
        pairs_path=pairs_path,
        output=output,
        side=side,
        end=end,
        shift=shift,
        window_size=window_size,
        chunksize=chunksize,
        threads=threads,
        chromsizes=chromsizes,
        output_bigwig=output_bigwig,
        **kwargs,
    )


if __name__ == "__main__":
    coverage()


def coverage_py(
    pairs_path,
    output,
    side,
    end,
    shift,
    window_size,
    chunksize,
    threads,
    chromsizes,
    output_bigwig,
    **kwargs,
):
    """Generate coverage from pairs file.

    PAIRS_PATH : input .pairs/.pairsam file. If the path ends with .gz or .lz4, the
    input is decompressed by bgzip/lz4c.
    By default, the input is read from stdin.
    """

    # Handle chromsizes - either from file or extract from pairs header
    # Read pairs and get chromsizes from header
    pairs_chunks, chromsizes_dict = read_pairs(
        pairs_path if pairs_path else "-", chunksize, threads
    )
    if chromsizes is not None:
        chromsizes_dict = bioframe.read_chromsizes(
            chromsizes, filter_chroms=False, natsort=False
        )

    # Calculate coverage using the API
    coverage_df = calculate_coverage(
        pairs_chunks=pairs_chunks,
        chromsizes=chromsizes_dict,
        side=side,
        end=end,
        shift=shift,
        window_size=window_size,
        threads=threads,
    )

    # Save the results
    coverage_df.to_csv(output, sep="\t", index=False, header=False)

    if output_bigwig is not None:
        bioframe.to_bigwig(
            coverage_df[["chrom", "start", "end", "count"]],
            chromsizes_dict,
            output_bigwig,
        )
