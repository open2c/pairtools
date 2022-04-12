#!/usr/bin/env python
import sys
import glob
import math
import subprocess
import click

from ..lib import fileio, pairsam_format, headerops
from . import cli, common_io_options

UTIL_NAME = "pairtools_merge"


@cli.command()
@click.argument(
    "pairs_path",
    nargs=-1,
    type=str,
)
@click.option(
    "-o",
    "--output",
    type=str,
    default="",
    help="output file."
    " If the path ends with .gz/.lz4, the output is compressed by bgzip/lz4c."
    " By default, the output is printed into stdout.",
)
@click.option(
    "--max-nmerge",
    type=int,
    default=8,
    show_default=True,
    help="The maximal number of inputs merged at once. For more, store "
    "merged intermediates in temporary files.",
)
@click.option(
    "--tmpdir",
    type=str,
    default="",
    help="Custom temporary folder for merged intermediates.",
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
    default="",
    show_default=True,
    help="A binary to compress temporary merged chunks. "
    "Must decompress input when the flag -d is provided. "
    "Suggested alternatives: lz4c, gzip, lzop, snzip. "
    "NOTE: fails silently if the command syntax is wrong. ",
)
@click.option(
    "--nproc",
    type=int,
    default=8,
    help="Number of threads for merging.",
    show_default=True,
)
@click.option(
    "--nproc-in",
    type=int,
    default=1,
    show_default=True,
    help="Number of processes used by the auto-guessed input decompressing command.",
)
@click.option(
    "--nproc-out",
    type=int,
    default=8,
    show_default=True,
    help="Number of processes used by the auto-guessed output compressing command.",
)
@click.option(
    "--cmd-in",
    type=str,
    default=None,
    help="A command to decompress the input. "
    "If provided, fully overrides the auto-guessed command. "
    "Does not work with stdin. "
    "Must read input from stdin and print output into stdout. "
    "EXAMPLE: pbgzip -dc -n 3",
)
@click.option(
    "--cmd-out",
    type=str,
    default=None,
    help="A command to compress the output. "
    "If provided, fully overrides the auto-guessed command. "
    "Does not work with stdout. "
    "Must read input from stdin and print output into stdout. "
    "EXAMPLE: pbgzip -c -n 8",
)
@click.option(
    "--keep-first-header/--no-keep-first-header",
    default=False,
    show_default=True,
    help="Keep the first header or merge the headers together. Default: merge headers.",
)
@click.option(
    "--concatenate/--no-concatenate",
    default=False,
    show_default=True,
    help="Simple concatenate instead of merging sorted files.",
)
# Using custom IO options


def merge(
    pairs_path, output, max_nmerge, tmpdir, memory, compress_program, nproc, **kwargs
):
    """Merge .pairs/.pairsam files.
    By default, assumes that the files are sorted and maintains the sorting.

    Merge triu-flipped sorted pairs/pairsam files. If present, the @SQ records
    of the SAM header must be identical; the sorting order of
    these lines is taken from the first file in the list.
    The ID fields of the @PG records of the SAM header are modified with a
    numeric suffix to produce unique records.
    The other unique SAM and non-SAM header lines are copied into the output header.

    PAIRS_PATH : upper-triangular flipped sorted .pairs/.pairsam files to merge
    or a group/groups of .pairs/.pairsam files specified by a wildcard. For
    paths ending in .gz/.lz4, the files are decompressed by bgzip/lz4c.

    """
    merge_py(
        pairs_path,
        output,
        max_nmerge,
        tmpdir,
        memory,
        compress_program,
        nproc,
        **kwargs,
    )


def merge_py(
    pairs_path, output, max_nmerge, tmpdir, memory, compress_program, nproc, **kwargs
):
    paths = sum([glob.glob(mask) for mask in pairs_path], [])

    if len(paths) == 0:
        raise ValueError(f"No input paths: {pairs_path}")

    outstream = fileio.auto_open(
        output,
        mode="w",
        nproc=kwargs.get("nproc_out"),
        command=kwargs.get("cmd_out", None),
    )

    # if there is only one input, bypass merging and do not modify the header
    if len(paths) == 1:
        instream = fileio.auto_open(
            paths[0],
            mode="r",
            nproc=kwargs.get("nproc_in"),
            command=kwargs.get("cmd_in", None),
        )
        for line in instream:
            outstream.write(line)
        if outstream != sys.stdout:
            outstream.close()

        return

    headers = []
    for path in paths:
        f = fileio.auto_open(
            path,
            mode="r",
            nproc=kwargs.get("nproc_in"),
            command=kwargs.get("cmd_in", None),
        )
        h, _ = headerops.get_header(f)
        headers.append(h)
        f.close()
        # Skip other headers if keep_first_header is True (False by default):
        if kwargs.get("keep_first_header", False):
            break

    if not headerops.all_same_columns(headers):
        raise ValueError("Input pairs cannot contain different columns")

    merged_header = headerops.merge_headers(headers)
    merged_header = headerops.append_new_pg(merged_header, ID=UTIL_NAME, PN=UTIL_NAME)

    outstream.writelines((l + "\n" for l in merged_header))
    outstream.flush()

    # If concatenation requested instead of merging sorted input:
    if kwargs.get("concatenate", False):
        command = r"""
                    /bin/bash -c 'export LC_COLLATE=C; export LANG=C; cat """
    # Full merge that keeps the ordered input:
    else:
        command = r"""
            /bin/bash -c 'export LC_COLLATE=C; export LANG=C; sort
            -k {0},{0} -k {1},{1} -k {2},{2}n -k {3},{3}n -k {4},{4} 
            --merge  
            --field-separator=$'\''{5}'\''
            {6}
            {7}
            {8}
            -S {9}
            {10}
            """.replace(
            "\n", " "
        ).format(
            pairsam_format.COL_C1 + 1,
            pairsam_format.COL_C2 + 1,
            pairsam_format.COL_P1 + 1,
            pairsam_format.COL_P2 + 1,
            pairsam_format.COL_PTYPE + 1,
            pairsam_format.PAIRSAM_SEP_ESCAPE,
            " --parallel={} ".format(nproc) if nproc > 1 else " ",
            " --batch-size={} ".format(max_nmerge) if max_nmerge else " ",
            " --temporary-directory={} ".format(tmpdir) if tmpdir else " ",
            memory,
            (
                " --compress-program={} ".format(compress_program)
                if compress_program
                else " "
            ),
        )
    for path in paths:
        if kwargs.get("cmd_in", None):
            command += r""" <(cat {} | {} | sed -n -e '\''/^[^#]/,$p'\'')""".format(
                path, kwargs["cmd_in"]
            )
        elif path.endswith(".gz"):
            command += (
                r""" <(bgzip -dc -@ {} {} | sed -n -e '\''/^[^#]/,$p'\'')""".format(
                    kwargs["nproc_in"], path
                )
            )
        elif path.endswith(".lz4"):
            command += r""" <(lz4c -dc {} | sed -n -e '\''/^[^#]/,$p'\'')""".format(
                path
            )
        else:
            command += r""" <(sed -n -e '\''/^[^#]/,$p'\'' {})""".format(path)
    command += "'"

    subprocess.check_call(command, shell=True, stdout=outstream)

    if outstream != sys.stdout:
        outstream.close()


if __name__ == "__main__":
    merge()
