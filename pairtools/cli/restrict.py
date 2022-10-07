#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import click
import warnings

import numpy as np

from ..lib import fileio, pairsam_format, headerops
from . import cli, common_io_options

from ..lib.restrict import find_rfrag

UTIL_NAME = "pairtools_restrict"


@cli.command()
@click.argument("pairs_path", type=str, required=False)
@click.option(
    "-f",
    "--frags",
    type=str,
    required=True,
    help="a tab-separated BED file with the positions of restriction fragments "
    "(chrom, start, end). Can be generated using cooler digest.",
)
@click.option(
    "-o",
    "--output",
    type=str,
    default="",
    help="output .pairs/.pairsam file."
    " If the path ends with .gz/.lz4, the output is compressed by bgzip/lz4c."
    " By default, the output is printed into stdout.",
)
@common_io_options
def restrict(pairs_path, frags, output, **kwargs):
    """Assign restriction fragments to pairs.

    Identify the restriction fragments that got ligated into a Hi-C molecule.

    Note: rfrags are 0-indexed

    PAIRS_PATH : input .pairs/.pairsam file. If the path ends with .gz/.lz4, the
    input is decompressed by bgzip/lz4c. By default, the input is read from stdin.
    """
    restrict_py(pairs_path, frags, output, **kwargs)


def restrict_py(pairs_path, frags, output, **kwargs):
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
    header = headerops.append_columns(
        header,
        [
            "rfrag1",
            "rfrag_start1",
            "rfrag_end1",
            "rfrag2",
            "rfrag_start2",
            "rfrag_end2",
        ],
    )

    outstream.writelines((l + "\n" for l in header))
    rfrags = np.genfromtxt(
        frags,
        delimiter="\t",
        comments="#",
        dtype=None,
        encoding="ascii",
        names=["chrom", "start", "end"],
    )

    rfrags.sort(order=["chrom", "start", "end"])
    chrom_borders = np.r_[
        0, 1 + np.where(rfrags["chrom"][:-1] != rfrags["chrom"][1:])[0], rfrags.shape[0]
    ]
    rfrags = {
        rfrags["chrom"][i]: np.concatenate([[0], rfrags["end"][i:j] + 1])
        for i, j in zip(chrom_borders[:-1], chrom_borders[1:])
    }

    for line in body_stream:
        cols = line.rstrip('\n').split(pairsam_format.PAIRSAM_SEP)
        chrom1, pos1 = cols[pairsam_format.COL_C1], int(cols[pairsam_format.COL_P1])
        rfrag_idx1, rfrag_start1, rfrag_end1 = find_rfrag(rfrags, chrom1, pos1)
        chrom2, pos2 = cols[pairsam_format.COL_C2], int(cols[pairsam_format.COL_P2])
        rfrag_idx2, rfrag_start2, rfrag_end2 = find_rfrag(rfrags, chrom2, pos2)
        cols += [str(rfrag_idx1), str(rfrag_start1), str(rfrag_end1)]
        cols += [str(rfrag_idx2), str(rfrag_start2), str(rfrag_end2)]
        outstream.write(pairsam_format.PAIRSAM_SEP.join(cols))
        outstream.write("\n")

    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()


if __name__ == "__main__":
    restrict()
