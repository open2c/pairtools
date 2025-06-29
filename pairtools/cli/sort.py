#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import sys
import click
import subprocess
import shutil
import warnings
from .._logging import get_logger

from ..lib import fileio, pairsam_format, headerops
from . import cli, common_io_options

logger = get_logger()
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
    "--c1",
    type=str,
    default=pairsam_format.COLUMNS_PAIRS[1],
    help=f"Chrom 1 column; default {pairsam_format.COLUMNS_PAIRS[1]}"
    "[input format option]",
)
@click.option(
    "--c2",
    type=str,
    default=pairsam_format.COLUMNS_PAIRS[3],
    help=f"Chrom 2 column; default {pairsam_format.COLUMNS_PAIRS[3]}"
    "[input format option]",
)
@click.option(
    "--p1",
    type=str,
    default=pairsam_format.COLUMNS_PAIRS[2],
    help=f"Position 1 column; default {pairsam_format.COLUMNS_PAIRS[2]}"
    "[input format option]",
)
@click.option(
    "--p2",
    type=str,
    default=pairsam_format.COLUMNS_PAIRS[4],
    help=f"Position 2 column; default {pairsam_format.COLUMNS_PAIRS[4]}"
    "[input format option]",
)
@click.option(
    "--pt",
    type=str,
    default=pairsam_format.COLUMNS_PAIRS[7],
    help=f"Pair type column; default {pairsam_format.COLUMNS_PAIRS[7]}"
    "[input format option]",
    required=False,
)
@click.option(
    "--extra-col",
    nargs=1,
    type=str,
    multiple=True,
    help="Extra column (name or numerical index) that is also used for sorting."
    "The option can be provided multiple times."
    'Example: --extra-col "phase1" --extra-col "phase2". [output format option]',
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
def sort(
    pairs_path,
    output,
    c1,
    c2,
    p1,
    p2,
    pt,
    extra_col,
    nproc,
    tmpdir,
    memory,
    compress_program,
    **kwargs,
):
    """Sort a .pairs/.pairsam file.

    Sort pairs in the lexicographic order along chrom1 and chrom2, in the
    numeric order along pos1 and pos2 and in the lexicographic order along
    pair_type.

    PAIRS_PATH : input .pairs/.pairsam file. If the path ends with .gz or .lz4, the
    input is decompressed by bgzip or lz4c, correspondingly. By default, the
    input is read as text from stdin.
    """
    sort_py(
        pairs_path,
        output,
        c1,
        c2,
        p1,
        p2,
        pt,
        extra_col,
        nproc,
        tmpdir,
        memory,
        compress_program,
        **kwargs,
    )


def sort_py(
    pairs_path,
    output,
    c1,
    c2,
    p1,
    p2,
    pt,
    extra_col,
    nproc,
    tmpdir,
    memory,
    compress_program,
    **kwargs,
):
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

    column_names = headerops.extract_column_names(header)
    column_names = headerops.canonicalize_columns(column_names)

    # Get column indices with fallbacks
    try:
        col1 = headerops.get_column_index(column_names, c1)
        col2 = headerops.get_column_index(column_names, c2)
        colp1 = headerops.get_column_index(column_names, p1)
        colp2 = headerops.get_column_index(column_names, p2)
        
        # Make pair_type optional
        try:
            colpt = headerops.get_column_index(column_names, pt) if pt else None
        except ValueError:
            colpt = None
            
        extra_cols = []
        for col in extra_col:
            try:
                extra_cols.append(headerops.get_column_index(column_names, col))
            except ValueError:
                logger.warning(f"Extra column {col} not found in header, skipping")
                continue
    except ValueError as e:
        raise ValueError(f"Column error: {str(e)}") from e

    # Generate sort command columns
    cols = []
    for i, col in enumerate([col1, colp1, col2, colp2, colpt] + extra_cols):
        if col is None:
            continue  # Skip optional columns that weren't found
        dtype = pairsam_format.DTYPES_PAIRSAM.get(column_names[col], str)
        cols.append(
            f"-k {col+1},{col+1}{'n' if issubclass(dtype, int) else ''}"
        )
    cols = " ".join(cols)
    
    command = rf"""
        /bin/bash -c 'export LC_COLLATE=C; export LANG=C; sort 
        {cols}
        --stable
        --field-separator=$'\''{pairsam_format.PAIRSAM_SEP_ESCAPE}'\''
        --parallel={nproc}
        {f'--temporary-directory={tmpdir}' if tmpdir else ''}
        -S {memory}
        {f'--compress-program={compress_program}' if compress_program else ''}'
        """.replace("\n", " ")
        
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
