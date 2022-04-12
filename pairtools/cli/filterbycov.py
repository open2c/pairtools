#!/usr/bin/env python
# -*- coding: utf-8  -*-
import sys
import ast
import warnings
import pathlib

import click

from ..lib import fileio, pairsam_format, headerops, dedup
from . import cli, common_io_options

from ..lib.filterbycov import streaming_filterbycov
from ..lib.stats import PairCounter


UTIL_NAME = "pairtools_filterbycov"

######################################
## TODO: - output stats after filtering
## edit/update mark as dup to mark as multi
###################################


@cli.command()
@click.argument("pairs_path", type=str, required=False)
@click.option(
    "-o",
    "--output",
    type=str,
    default="",
    help="output file for pairs from low coverage regions."
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4c-compressed."
    " By default, the output is printed into stdout.",
)
@click.option(
    "--output-highcov",
    type=str,
    default="",
    help="output file for pairs from high coverage regions."
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4c-compressed."
    " If the path is the same as in --output or -, output duplicates together "
    " with deduped pairs. By default, duplicates are dropped.",
)
@click.option(
    "--output-unmapped",
    type=str,
    default="",
    help="output file for unmapped pairs. "
    "If the path ends with .gz or .lz4, the output is bgzip-/lz4c-compressed. "
    "If the path is the same as in --output or -, output unmapped pairs together "
    "with deduped pairs. If the path is the same as --output-highcov, "
    "output unmapped reads together. By default, unmapped pairs are dropped.",
)
@click.option(
    "--output-stats",
    type=str,
    default="",
    help="output file for statistics of multiple interactors. "
    " If file exists, it will be open in the append mode."
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4c-compressed."
    " By default, statistics are not printed.",
)
@click.option(
    "--max-cov", type=int, default=8, help="The maximum allowed coverage per region."
)
@click.option(
    "--max-dist",
    type=int,
    default=500,
    help="The resolution for calculating coverage. For each pair, the local "
    "coverage around each end is calculated as (1 + the number of neighbouring "
    "pairs within +/- max_dist bp) ",
)
@click.option(
    "--method",
    type=click.Choice(["max", "sum"]),
    default="max",
    help="calculate the number of neighbouring pairs as either the sum or the max"
    " of the number of neighbours on the two sides",
    show_default=True,
)
@click.option(
    "--sep",
    type=str,
    default=pairsam_format.PAIRSAM_SEP_ESCAPE,
    help=r"Separator (\t, \v, etc. characters are " "supported, pass them in quotes) ",
)
@click.option(
    "--comment-char", type=str, default="#", help="The first character of comment lines"
)
@click.option(
    "--send-header-to",
    type=click.Choice(["lowcov", "highcov", "both", "none"]),
    default="both",
    help="Which of the outputs should receive header and comment lines",
)
@click.option(
    "--c1",
    type=int,
    default=pairsam_format.COL_C1,
    help="Chrom 1 column; default {}".format(pairsam_format.COL_C1),
)
@click.option(
    "--c2",
    type=int,
    default=pairsam_format.COL_C2,
    help="Chrom 2 column; default {}".format(pairsam_format.COL_C2),
)
@click.option(
    "--p1",
    type=int,
    default=pairsam_format.COL_P1,
    help="Position 1 column; default {}".format(pairsam_format.COL_P1),
)
@click.option(
    "--p2",
    type=int,
    default=pairsam_format.COL_P2,
    help="Position 2 column; default {}".format(pairsam_format.COL_P2),
)
@click.option(
    "--s1",
    type=int,
    default=pairsam_format.COL_S1,
    help="Strand 1 column; default {}".format(pairsam_format.COL_S1),
)
@click.option(
    "--s2",
    type=int,
    default=pairsam_format.COL_S2,
    help="Strand 2 column; default {}".format(pairsam_format.COL_S2),
)
@click.option(
    "--unmapped-chrom",
    type=str,
    default=pairsam_format.UNMAPPED_CHROM,
    help="Placeholder for a chromosome on an unmapped side; default {}".format(
        pairsam_format.UNMAPPED_CHROM
    ),
)
@click.option(
    "--mark-multi",
    is_flag=True,
    help='If specified, duplicate pairs are marked as FF in "pair_type" and '
    "as a duplicate in the sam entries.",
)
@common_io_options
def filterbycov(
    pairs_path,
    output,
    output_highcov,
    output_unmapped,
    output_stats,
    max_dist,
    max_cov,
    method,
    sep,
    comment_char,
    send_header_to,
    c1,
    c2,
    p1,
    p2,
    s1,
    s2,
    unmapped_chrom,
    mark_multi,
    **kwargs
):
    """Remove pairs from regions of high coverage.

    Find and remove pairs with >(MAX_COV-1) neighbouring pairs
    within a +/- MAX_DIST bp window around either side. Useful for single-cell
    Hi-C experiments, where coverage is naturally limited by the chromosome
    copy number.

    PAIRS_PATH : input triu-flipped sorted .pairs or .pairsam file.  If the
    path ends with .gz/.lz4, the input is decompressed by bgzip/lz4c.
    By default, the input is read from stdin.
    """
    filterbycov_py(
        pairs_path,
        output,
        output_highcov,
        output_unmapped,
        output_stats,
        max_dist,
        max_cov,
        method,
        sep,
        comment_char,
        send_header_to,
        c1,
        c2,
        p1,
        p2,
        s1,
        s2,
        unmapped_chrom,
        mark_multi,
        **kwargs
    )


def filterbycov_py(
    pairs_path,
    output,
    output_highcov,
    output_unmapped,
    output_stats,
    max_dist,
    max_cov,
    method,
    sep,
    comment_char,
    send_header_to,
    c1,
    c2,
    p1,
    p2,
    s1,
    s2,
    unmapped_chrom,
    mark_multi,
    **kwargs
):

    ## Prepare input, output streams based on selected outputs
    ## Default ouput stream is low-frequency interactors
    sep = ast.literal_eval('"""' + sep + '"""')
    send_header_to_lowcov = send_header_to in ["both", "lowcov"]
    send_header_to_highcov = send_header_to in ["both", "highcov"]

    instream = (
        fileio.auto_open(
            pairs_path,
            mode="r",
            nproc=kwargs.get("nproc_in"),
            command=kwargs.get("cmd_in", None),
        )
        if pairs_path
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
    out_stats_stream = (
        fileio.auto_open(
            output_stats,
            mode="w",
            nproc=kwargs.get("nproc_out"),
            command=kwargs.get("cmd_out", None),
        )
        if output_stats
        else None
    )

    # generate empty PairCounter if stats output is requested:
    out_stat = PairCounter() if output_stats else None

    # output the high-frequency interacting pairs
    if not output_highcov:
        outstream_high = None
    elif output_highcov == "-" or (
        pathlib.Path(output_highcov).absolute() == pathlib.Path(output).absolute()
    ):
        outstream_high = outstream
    else:
        outstream_high = fileio.auto_open(
            output_highcov,
            mode="w",
            nproc=kwargs.get("nproc_out"),
            command=kwargs.get("cmd_out", None),
        )

    # output unmapped pairs
    if not output_unmapped:
        outstream_unmapped = None
    elif output_unmapped == "-" or (
        pathlib.Path(output_unmapped).absolute() == pathlib.Path(output).absolute()
    ):
        outstream_unmapped = outstream
    elif (
        pathlib.Path(output_unmapped).absolute()
        == pathlib.Path(output_highcov).absolute()
    ):
        outstream_unmapped = outstream_high
    else:
        outstream_unmapped = fileio.auto_open(
            output_unmapped,
            mode="w",
            nproc=kwargs.get("nproc_out"),
            command=kwargs.get("cmd_out", None),
        )

    # prepare file headers
    header, body_stream = headerops.get_header(instream)
    header = headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)

    # header for low-frequency interactors
    if send_header_to_lowcov:
        outstream.writelines((l + "\n" for l in header))

    # header for high-frequency interactors
    if send_header_to_highcov and outstream_high and (outstream_high != outstream):
        outstream_high.writelines((l + "\n" for l in header))

    # header for unmapped pairs
    if (
        outstream_unmapped
        and (outstream_unmapped != outstream)
        and (outstream_unmapped != outstream_high)
    ):
        outstream_unmapped.writelines((l + "\n" for l in header))

    # perform filtering of pairs based on low/high-frequency of interaction
    streaming_filterbycov(
        method,
        max_dist,
        max_cov,
        sep,
        c1,
        c2,
        p1,
        p2,
        s1,
        s2,
        unmapped_chrom,
        body_stream,
        outstream,
        outstream_high,
        outstream_unmapped,
        out_stat,
        mark_multi,
    )

    ## FINISHED!
    # save statistics to a file if it was requested: TO BE TESTED
    if out_stat:
        out_stat.save(out_stats_stream)

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()

    if outstream_high and (outstream_high != outstream):
        outstream_high.close()

    if (
        outstream_unmapped
        and (outstream_unmapped != outstream)
        and (outstream_unmapped != outstream_high)
    ):
        outstream_unmapped.close()

    if out_stats_stream:
        out_stats_stream.close()


if __name__ == "__main__":
    filterbycov()
