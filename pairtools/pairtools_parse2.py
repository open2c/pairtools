#!/usr/bin/env python
# -*- coding: utf-8 -*-
from collections import OrderedDict
import subprocess
import fileinput
import itertools
import click
import pipes
import sys
import os
import io

from . import _fileio, _pairsam_format, _parse2, _headerops, cli, common_io_options
from .pairtools_stats import PairCounter

UTIL_NAME = "pairtools_parse2"

EXTRA_COLUMNS = [
    "mapq",
    "pos5",
    "pos3",
    "cigar",
    "read_len",
    "matched_bp",
    "algn_ref_span",
    "algn_read_span",
    "dist_to_5",
    "dist_to_3",
    "seq",
]


@cli.command()
@click.argument("sam_path", type=str, required=False)

# Parsing options:
@click.option(
    "-c",
    "--chroms-path",
    type=str,
    required=True,
    help="Chromosome order used to flip interchromosomal mates: "
    "path to a chromosomes file (e.g. UCSC chrom.sizes or similar) whose "
    "first column lists scaffold names. Any scaffolds not listed will be "
    "ordered lexicographically following the names provided.",
)
@click.option(
    "--assembly",
    type=str,
    help="Name of genome assembly (e.g. hg19, mm10) to store in the pairs header.",
)
@click.option(
    "--min-mapq",
    type=int,
    default=1,
    show_default=True,
    help="The minimal MAPQ score to consider a read as uniquely mapped",
)
@click.option(
    "--max-inter-align-gap",
    type=int,
    default=20,
    show_default=True,
    help="read segments that are not covered by any alignment and"
    ' longer than the specified value are treated as "null" alignments.'
    " These null alignments convert otherwise linear alignments into walks,"
    " and affect how they get reported as a Hi-C pair.",
)
@click.option(
    "--max-fragment-size",
    type=int,
    default=500,
    show_default=True,
    help="Largest fragment size for the detection of overlapping "
    "alignments at the ends of forward and reverse reads. "
    "Not used in --single-end mode. ",
)
@click.option(
    "--single-end", is_flag=True, help="If specified, the input is single-end."
)

# Reporting options:
@click.option(
    "-o",
    "--output-file",
    type=str,
    default="",
    help="output file. "
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4-compressed."
    "By default, the output is printed into stdout. ",
)
@click.option(
    "--coordinate-system",
    type=click.Choice(["read", "walk", "pair"]),
    default="read",
    help="coordinate system for reporting the walk. "
    ' "read" - orient each pair as it appeared on a read, starting from 5\'-end of forward then reverse read. '
    ' "walk" - orient each pair as it appeared sequentially in the reconstructed walk. '
    ' "pair" - re-orient each pair as if it was sequenced independently by Hi-C. ',
    show_default=True,
)
@click.option(
    "--no-flip",
    is_flag=True,
    help="If specified, do not flip pairs in genomic order and instead preserve "
    "the order in which they were sequenced.",
)
@click.option(
    "--drop-readid",
    is_flag=True,
    help="If specified, do not add read ids to the output",
)
@click.option(
    "--readid-transform",
    type=str,
    default=None,
    help="A Python expression to modify read IDs. Useful when read IDs differ "
    "between the two reads of a pair. Must be a valid Python expression that "
    "uses variables called readID and/or i (the 0-based index of the read pair "
    "in the bam file) and returns a new value, e.g. \"readID[:-2]+'_'+str(i)\". "
    "Make sure that transformed readIDs remain unique!",
    show_default=True,
)
@click.option(
    "--drop-seq",
    is_flag=True,
    help="If specified, remove sequences and PHREDs from the sam fields",
)
@click.option(
    "--drop-sam", is_flag=True, help="If specified, do not add sams to the output"
)
@click.option(
    "--add-junction-index",
    is_flag=True,
    help="If specified, parse2 will report junction index for each pair in the walk",
)
@click.option(
    "--add-columns",
    type=click.STRING,
    default="",
    help="Report extra columns describing alignments "
    "Possible values (can take multiple values as a comma-separated "
    "list): a SAM tag (any pair of uppercase letters) or {}.".format(
        ", ".join(EXTRA_COLUMNS)
    ),
)
@click.option(
    "--output-stats",
    type=str,
    default="",
    help="output file for various statistics of pairs file. "
    " By default, statistics is not generated.",
)
@common_io_options
def parse2(
    sam_path,
    chroms_path,
    output_file,
    assembly,
    min_mapq,
    drop_readid,
    drop_seq,
    drop_sam,
    add_junction_index,
    add_columns,
    output_stats,
    coordinate_system,
    **kwargs
):
    """Find ligation junctions in .sam, make .pairs.
    SAM_PATH : an input .sam/.bam file with paired-end sequence alignments of
    Hi-C molecules. If the path ends with .bam, the input is decompressed from
    bam with samtools. By default, the input is read from stdin.
    """
    parse2_py(
        sam_path,
        chroms_path,
        output_file,
        assembly,
        min_mapq,
        drop_readid,
        drop_seq,
        drop_sam,
        add_junction_index,
        add_columns,
        output_stats,
        coordinate_system,
        **kwargs
    )


def parse2_py(
    sam_path,
    chroms_path,
    output_file,
    assembly,
    min_mapq,
    drop_readid,
    drop_seq,
    drop_sam,
    add_junction_index,
    add_columns,
    output_stats,
    coordinate_system,
    **kwargs
):
    instream = (
        _fileio.auto_open(
            sam_path,
            mode="r",
            nproc=kwargs.get("nproc_in"),
            command=kwargs.get("cmd_in", None),
        )
        if sam_path
        else sys.stdin
    )
    outstream = (
        _fileio.auto_open(
            output_file,
            mode="w",
            nproc=kwargs.get("nproc_out"),
            command=kwargs.get("cmd_out", None),
        )
        if output_file
        else sys.stdout
    )
    out_stats_stream = (
        _fileio.auto_open(
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

    samheader, body_stream = _headerops.get_header(instream, comment_char="@")

    if not samheader:
        raise ValueError(
            "The input sam is missing a header! If reading a bam file, please use `samtools view -h` to include the header."
        )

    sam_chromsizes = _headerops.get_chromsizes_from_sam_header(samheader)
    chromosomes = _headerops.get_chrom_order(chroms_path, list(sam_chromsizes.keys()))

    add_columns = [col for col in add_columns.split(",") if col]
    for col in add_columns:
        if not ((col in EXTRA_COLUMNS) or (len(col) == 2 and col.isupper())):
            raise Exception("{} is not a valid extra column".format(col))

    columns = _pairsam_format.COLUMNS + (
        [c + side for c in add_columns for side in ["1", "2"]]
    )

    if drop_sam:
        columns.pop(columns.index("sam1"))
        columns.pop(columns.index("sam2"))

    if not add_junction_index:
        columns.pop(columns.index("junction_index"))

    header = _headerops.make_standard_pairsheader(
        assembly=assembly,
        chromsizes=[(chrom, sam_chromsizes[chrom]) for chrom in chromosomes],
        columns=columns,
        shape="whole matrix" if coordinate_system != "pair" else "upper triangle",
    )

    header = _headerops.insert_samheader(header, samheader)
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    outstream.writelines((l + "\n" for l in header))

    _parse2.streaming_classify(
        body_stream,
        outstream,
        chromosomes,
        min_mapq,
        drop_readid,
        drop_seq,
        drop_sam,
        add_junction_index,
        add_columns,
        out_stat,
        coordinate_system,
        **kwargs
    )

    # save statistics to a file if it was requested:
    if out_stat:
        out_stat.save(out_stats_stream)

    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()
    if out_stats_stream:
        out_stats_stream.close()
