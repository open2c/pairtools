#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import sys
import click
import pandas as pd

from ..lib import fileio
from . import cli, common_io_options

from ..lib.scaling import compute_scaling


UTIL_NAME = "pairtools_scaling"


@cli.command()
@click.argument("input_path", type=str, nargs=-1, required=False)
@click.option(
    "-o", "--output", type=str, default="", help="output .tsv file with summary."
)
@click.option(
    "--view",
    "--regions",
    help="Path to a BED file which defines which regions (viewframe) of the chromosomes to use. "
    "By default, this is parsed from .pairs header. ",
    type=str,
    required=False,
    default=None,
)
@click.option(
    "--chunksize",
    type=int,
    default=100_000,
    show_default=True,
    required=False,
    help="Number of pairs in each chunk. Reduce for lower memory footprint.",
)
@click.option(
    "--dist-range",
    type=click.Tuple([int, int]),
    default=(1, 1_000_000_000),
    show_default=True,
    required=False,
    help="Distance range. ",
)
@click.option(
    "--n-dist-bins-decade",
    type=int,
    default=8,
    show_default=True,
    required=False,
    help="Number of bins to split the distance range in log10-space, specified per a factor of 10 difference.",
)
@common_io_options
def scaling(
    input_path, output, view, chunksize, dist_range, n_dist_bins_decade, **kwargs
):
    """Calculate pairs scalings.

    INPUT_PATH : by default, a .pairs/.pairsam file to calculate statistics.
    If not provided, the input is read from stdin.

    The files with paths ending with .gz/.lz4 are decompressed by bgzip/lz4c.

    Output is .tsv file with scaling stats (both cis scalings and trans levels).
    """
    scaling_py(
        input_path, output, view, chunksize, dist_range, n_dist_bins_decade, **kwargs
    )


def scaling_py(
    input_path, output, view, chunksize, dist_range, n_dist_bins_decade, **kwargs
):

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

    if view is not None:
        view = pd.read_table(view)

    # Pass the header to the instream so that it can parse the header automatically
    cis_scalings, trans_levels = compute_scaling(
        instream,
        regions=view,
        chromsizes=None,
        dist_range=dist_range,
        n_dist_bins_decade=n_dist_bins_decade,
        chunksize=chunksize,
    )
    summary_stats = pd.concat([cis_scalings, trans_levels])

    # save statistics to the file
    summary_stats.to_csv(outstream, sep="\t", index=False)

    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()


if __name__ == "__main__":
    scaling()
