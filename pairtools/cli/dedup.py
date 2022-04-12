#!/usr/bin/env python
# -*- coding: utf-8  -*-

import sys
import ast
import pathlib

# from distutils.log import warn
# import warnings

from ..lib import fileio, pairsam_format, headerops
from . import cli, common_io_options

import click

from ..lib.dedup import streaming_dedup, streaming_dedup_cython
from ..lib.stats import PairCounter

UTIL_NAME = "pairtools_dedup"


@cli.command()
@click.argument("pairs_path", type=str, required=False)

### Output files:
@click.option(
    "-o",
    "--output",
    type=str,
    default="",
    help="output file for pairs after duplicate removal."
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4c-compressed."
    " By default, the output is printed into stdout.",
)
@click.option(
    "--output-dups",
    type=str,
    default="",
    help="output file for duplicated pairs. "
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
    "with deduped pairs. If the path is the same as --output-dups, output "
    "unmapped reads together with dups. By default, unmapped pairs are dropped.",
)
@click.option(
    "--output-stats",
    type=str,
    default="",
    help="output file for duplicate statistics."
    " If file exists, it will be open in the append mode."
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4c-compressed."
    " By default, statistics are not printed.",
)

### Set the dedup method:
@click.option(
    "--max-mismatch",
    type=int,
    default=3,
    show_default=True,
    help="Pairs with both sides mapped within this distance (bp) from each "
    "other are considered duplicates. ",
)
@click.option(
    "--method",
    type=click.Choice(["max", "sum"]),
    default="max",
    help="define the mismatch as either the max or the sum of the mismatches of"
    "the genomic locations of the both sides of the two compared molecules",
    show_default=True,
)
@click.option(
    "--backend",
    type=click.Choice(["scipy", "sklearn", "cython"]),
    default="scipy",
    help="What backend to use: scipy and sklearn are based on KD-trees,"
    " cython is online indexed list-based algorithm."
    " With cython backend, duplication is not transitive with non-zero max mismatch "
    " (e.g. pairs A and B are duplicates, and B and C are duplicates, then A and C are "
    " not necessary duplicates of each other), while with scipy and sklearn it's "
    " transitive (i.e. A and C are necessarily duplicates)."
    " Cython is the original version used in pairtools since its beginning."
    " It is available for backwards compatibility and to allow specification of the"
    " column order."
    " Now the default scipy backend is generally the fastest, and with chunksize below"
    " 1 mln has the lowest memory requirements."
    # " 'cython' is deprecated and provided for backwards compatibility",
)

### Scipy and sklearn-specific options:
@click.option(
    "--chunksize",
    type=int,
    default=100_000,
    show_default=True,
    help="Number of pairs in each chunk. Reduce for lower memory footprint."
    " Below 10,000 performance starts suffering significantly and the algorithm might"
    " miss a few duplicates with non-zero --max-mismatch."
    " Only works with '--backend scipy or sklearn'",
)
@click.option(
    "--carryover",
    type=int,
    default=100,
    show_default=True,
    help="Number of deduped pairs to carry over from previous chunk to the new chunk"
    " to avoid breaking duplicate clusters."
    " Only works with '--backend scipy or sklearn'",
)
@click.option(
    "-p",
    "--n-proc",
    type=int,
    default=1,
    help="Number of cores to use. Only applies with sklearn backend."
    "Still needs testing whether it is ever useful.",
)

### Output options:
@click.option(
    "--mark-dups",
    is_flag=True,
    help='If specified, duplicate pairs are marked as DD in "pair_type" and '
    "as a duplicate in the sam entries.",
)
@click.option(
    "--keep-parent-id",
    is_flag=True,
    help="If specified, duplicate pairs are marked with the readID of the retained"
    " deduped read in the 'parent_readID' field.",
)
@click.option(
    "--extra-col-pair",
    nargs=2,
    # type=click.Tuple([str, str]),
    multiple=True,
    help="Extra columns that also must match for two pairs to be marked as "
    "duplicates. Can be either provided as 0-based column indices or as column "
    'names (requires the "#columns" header field). The option can be provided '
    "multiple times if multiple column pairs must match. "
    'Example: --extra-col-pair "phase1" "phase2"',
)

### Input options:
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
    type=click.Choice(["dups", "dedup", "both", "none"]),
    default="both",
    help="Which of the outputs should receive header and comment lines",
)
@click.option(
    "--c1",
    type=int,
    default=pairsam_format.COL_C1,
    help=f"Chrom 1 column; default {pairsam_format.COL_C1}"
    " Only works with '--backend cython'",
)
@click.option(
    "--c2",
    type=int,
    default=pairsam_format.COL_C2,
    help=f"Chrom 2 column; default {pairsam_format.COL_C2}"
    " Only works with '--backend cython'",
)
@click.option(
    "--p1",
    type=int,
    default=pairsam_format.COL_P1,
    help=f"Position 1 column; default {pairsam_format.COL_P1}"
    " Only works with '--backend cython'",
)
@click.option(
    "--p2",
    type=int,
    default=pairsam_format.COL_P2,
    help=f"Position 2 column; default {pairsam_format.COL_P2}"
    " Only works with '--backend cython'",
)
@click.option(
    "--s1",
    type=int,
    default=pairsam_format.COL_S1,
    help=f"Strand 1 column; default {pairsam_format.COL_S1}"
    " Only works with '--backend cython'",
)
@click.option(
    "--s2",
    type=int,
    default=pairsam_format.COL_S2,
    help=f"Strand 2 column; default {pairsam_format.COL_S2}"
    " Only works with '--backend cython'",
)
@click.option(
    "--unmapped-chrom",
    type=str,
    default=pairsam_format.UNMAPPED_CHROM,
    help="Placeholder for a chromosome on an unmapped side; default {}".format(
        pairsam_format.UNMAPPED_CHROM
    ),
)
@common_io_options
def dedup(
    pairs_path,
    output,
    output_dups,
    output_unmapped,
    output_stats,
    chunksize,
    carryover,
    max_mismatch,
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
    mark_dups,
    extra_col_pair,
    keep_parent_id,
    backend,
    n_proc,
    **kwargs,
):
    """Find and remove PCR/optical duplicates.

    Find PCR duplicates in an upper-triangular flipped sorted pairs/pairsam
    file. Allow for a +/-N bp mismatch at each side of duplicated molecules.

    PAIRS_PATH : input triu-flipped sorted .pairs or .pairsam file.  If the
    path ends with .gz/.lz4, the input is decompressed by bgzip/lz4c.
    By default, the input is read from stdin.
    """

    dedup_py(
        pairs_path,
        output,
        output_dups,
        output_unmapped,
        output_stats,
        chunksize,
        carryover,
        max_mismatch,
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
        mark_dups,
        extra_col_pair,
        keep_parent_id,
        backend,
        n_proc,
        **kwargs,
    )


if __name__ == "__main__":
    dedup()


def dedup_py(
    pairs_path,
    output,
    output_dups,
    output_unmapped,
    output_stats,
    chunksize,
    carryover,
    max_mismatch,
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
    mark_dups,
    extra_col_pair,
    keep_parent_id,
    backend,
    n_proc,
    **kwargs,
):

    sep = ast.literal_eval('"""' + sep + '"""')
    send_header_to_dedup = send_header_to in ["both", "dedup"]
    send_header_to_dup = send_header_to in ["both", "dups"]

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

    if not output_dups:
        outstream_dups = None
    elif output_dups == "-" or (
        pathlib.Path(output_dups).absolute() == pathlib.Path(output).absolute()
    ):
        outstream_dups = outstream
    else:
        outstream_dups = fileio.auto_open(
            output_dups,
            mode="w",
            nproc=kwargs.get("nproc_out"),
            command=kwargs.get("cmd_out", None),
        )

    if not output_unmapped:
        outstream_unmapped = None
    elif output_unmapped == "-" or (
        pathlib.Path(output_unmapped).absolute() == pathlib.Path(output).absolute()
    ):
        outstream_unmapped = outstream
    elif (
        pathlib.Path(output_unmapped).absolute() == pathlib.Path(output_dups).absolute()
    ):
        outstream_unmapped = outstream_dups
    else:
        outstream_unmapped = fileio.auto_open(
            output_unmapped,
            mode="w",
            nproc=kwargs.get("nproc_out"),
            command=kwargs.get("cmd_out", None),
        )

    header, body_stream = headerops.get_header(instream)
    header = headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    if send_header_to_dedup:
        outstream.writelines((l + "\n" for l in header))
    if send_header_to_dup and outstream_dups and (outstream_dups != outstream):
        dups_header = header
        if keep_parent_id and len(dups_header) > 0:
            dups_header = headerops.append_columns(dups_header, ["parent_readID"])
        outstream_dups.writelines((l + "\n" for l in dups_header))
    if (
        outstream_unmapped
        and (outstream_unmapped != outstream)
        and (outstream_unmapped != outstream_dups)
    ):
        outstream_unmapped.writelines((l + "\n" for l in header))

    column_names = headerops.extract_column_names(header)
    extra_cols1 = []
    extra_cols2 = []
    if extra_col_pair is not None:
        for col1, col2 in extra_col_pair:
            extra_cols1.append(column_names[col1] if col1.isdigit() else col1)
            extra_cols2.append(column_names[col2] if col2.isdigit() else col2)

    if backend == "cython":
        # warnings.warn(
        #     "'cython' backend is deprecated and provided only"
        #     " for backwards compatibility",
        #     DeprecationWarning,
        # )
        extra_cols1 = [column_names.index(col) for col in extra_cols1]
        extra_cols2 = [column_names.index(col) for col in extra_cols2]
        streaming_dedup_cython(
            method,
            max_mismatch,
            sep,
            c1,
            c2,
            p1,
            p2,
            s1,
            s2,
            extra_cols1,
            extra_cols2,
            unmapped_chrom,
            body_stream,
            outstream,
            outstream_dups,
            outstream_unmapped,
            out_stat,
            mark_dups,
            keep_parent_id,
        )
    elif backend in ("scipy", "sklearn"):
        streaming_dedup(
            in_stream=instream,
            colnames=column_names,
            chunksize=chunksize,
            carryover=carryover,
            method=method,
            mark_dups=mark_dups,
            max_mismatch=max_mismatch,
            extra_col_pairs=list(extra_col_pair),
            keep_parent_id=keep_parent_id,
            unmapped_chrom=unmapped_chrom,
            comment_char=comment_char,
            outstream=outstream,
            outstream_dups=outstream_dups,
            outstream_unmapped=outstream_unmapped,
            out_stat=out_stat,
            backend=backend,
            n_proc=n_proc,
        )
    else:
        raise ValueError("Unknown backend")

    # save statistics to a file if it was requested:
    if out_stat:
        out_stat.save(out_stats_stream)

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()

    if outstream_dups and (outstream_dups != outstream):
        outstream_dups.close()

    if (
        outstream_unmapped
        and (outstream_unmapped != outstream)
        and (outstream_unmapped != outstream_dups)
    ):
        outstream_unmapped.close()

    if out_stats_stream:
        out_stats_stream.close()
