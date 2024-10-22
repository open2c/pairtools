# !/usr/bin/env python
# -*- coding: utf-8 -*-

import click
import sys

from ..lib import fileio, pairsam_format, headerops
from . import cli, common_io_options

from ..lib.stats import PairCounter
from ..lib.parse_pysam import AlignmentFilePairtoolized
from ..lib.parse import streaming_classify

UTIL_NAME = "pairtools_parse2"


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
    "-o",
    "--output",
    type=str,
    default="",
    help="output file with pairs. "
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4-compressed."
    "By default, the output is printed into stdout. ",
)
@click.option(
    "--report-position",
    type=click.Choice(["junction", "read", "walk", "outer"]),
    default="outer",
    help="""Reported position of alignments in pairs of complex walks (pos columns). 
    Each alignment in .bam/.sam Hi-C-like data has two ends, and you can report one or another depending of the position of alignment on a read or in a pair.  
    
    "junction" - inner ends of sequential alignments in each pair, aka ligation junctions,
    "read" - 5'-end of alignments relative to R1 or R2 read coordinate system (as in traditional Hi-C),
    "walk" - 5'-end of alignments relative to the whole walk coordinate system,
    "outer" - outer ends of sequential alignments in each pair (parse2 default). """,
)
@click.option(
    "--report-orientation",
    type=click.Choice(["pair", "read", "walk", "junction"]),
    default="pair",
    help="""Reported orientataion of pairs in complex walk (strand columns).
    Each alignment in .bam/.sam Hi-C-like data has orientation, and you can report it relative to the read, pair or whole walk coordinate system.
    
    "pair" - orientation as if each pair in complex walk was sequenced independently from the outer ends or molecule (as in traditional Hi-C, also complex walks default),
    "read" - orientation defined by the read (R1 or R2 read coordinate system),
    "walk" - orientation defined by the walk coordinate system,
    "junction" - reversed "pair" orientation, as if pair was sequenced in both directions starting from the junction""",
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
    help="The minimal MAPQ score to consider a read as uniquely mapped.",
)
@click.option(
    "--max-inter-align-gap",
    type=int,
    default=20,
    show_default=True,
    help="Read segments that are not covered by any alignment and"
    ' longer than the specified value are treated as "null" alignments.'
    " These null alignments convert otherwise linear alignments into walks,"
    " and affect how they get reported as a Hi-C pair.",
)
@click.option(
    "--max-insert-size",
    type=int,
    default=500,
    show_default=True,
    help="When searching for overlapping ends of left and right read (R1 and R2), this sets the minimal distance "
    "when two alignments on the same strand and chromosome are considered part of the same fragment (and thus reported as the same alignment "
    "and not a pair). For traditional Hi-C with long restriction fragments and shorter molecules after ligation+sonication, this "
    "can be the expected molecule size. For complex walks with short restriction fragments, this can be the expected restriction fragment "
    "size. Note that unsequenced insert is *terra incognita* and might contain unsequenced DNA (including ligations) in it. "
    "This parameter is ignored in --single-end mode. ",
)
@click.option(
    "--dedup-max-mismatch",
    type=int,
    default=3,
    show_default=True,
    help="Allowed mismatch between intramolecular alignments to detect readthrough duplicates. "
    "Pairs with both sides mapped within this distance (bp) from each "
    "other are considered duplicates. ",
)
@click.option(
    "--single-end",
    is_flag=True,
    help="If specified, the input is single-end. "
    "Never use this for paired-end data, because R1 read will be omitted. "
    "If single-end data is provided, but parameter is unset, the pairs will be "
    "generated, but may contain artificial UN pairs. ",
)
@click.option(
    "--expand/--no-expand",
    is_flag=True,
    help="If specified, perform combinatorial expansion on the pairs. "
    "Combinatorial expansion is a way to increase the number of contacts in you data, assuming that all DNA fragments in the same molecule (read) are in contact. "
    "Expanded pairs have modified pair type, 'E{separation}_{pair type}'",
)
@click.option(
    "--max-expansion-depth",
    type=int,
    default=None,
    show_default=True,
    help="Works in combination with --expand. "
    "Maximum number of segments separating pair. By default, expanding all possible pairs."
    "Setting the number will limit the expansion depth and enforce contacts from the same "
    "side of the read. ",
)
@click.option(
    "--add-pair-index",
    is_flag=True,
    help="If specified, parse2 will report pair index in the walk as additional columns (R1, R2, R1&R2 or R1-R2). "
    "See documentation: https://pairtools.readthedocs.io/en/latest/parsing.html#rescuing-complex-walks "
    "For combinatorial expanded pairs, two numbers will be reported: "
    "original pair index of the left and right segments. ",
)
@click.option(
    "--flip/--no-flip",
    is_flag=True,
    default=False,
    help="If specified, flip pairs in genomic order and instead preserve "
    "the order in which they were sequenced. Note that no flip is recommended for analysis of walks because it will "
    "override the order of alignments in pairs. Flip is required for appropriate deduplication of sorted pairs. "
    "Flip is not required for cooler cload, which runs flipping internally. ",
)
@click.option(
    "--add-columns",
    type=click.STRING,
    default="",
    help="Report extra columns describing alignments "
    "Possible values (can take multiple values as a comma-separated "
    "list): a SAM tag (any pair of uppercase letters) or {}.".format(
        ", ".join(pairsam_format.EXTRA_COLUMNS)
    ),
)
@click.option(
    "--drop-readid/--keep-readid",
    is_flag=True,
    default=False,
    help="If specified, do not add read ids to the output. By default, keep read ids. Useful for long walks analysis. ",
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
    "--drop-seq/--keep-seq",
    is_flag=True,
    default=False,
    help="Remove sequences and PHREDs from the sam fields by default. Kept otherwise. ",
)
@click.option(
    "--drop-sam/--keep-sam",
    is_flag=True,
    default=False,
    help="Do not add sams to the output by default. Kept otherwise. ",
)
@click.option(
    "--output-parsed-alignments",
    type=str,
    default="",
    help="output file with all parsed alignments (one alignment per line)."
    " Useful for debugging and analysis of walks."
    " If file exists, it will be open in the append mode."
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4-compressed."
    " By default, not used.",
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
    sam_path, chroms_path, output, output_parsed_alignments, output_stats, **kwargs
):
    """Extracts pairs from .sam/.bam data with complex walks, make .pairs.
    SAM_PATH : an input .sam/.bam file with paired-end or single-end sequence alignments of
    Hi-C (or Hi-C-like) molecules. If the path ends with .bam, the input is decompressed from
    bam with samtools. By default, the input is read from stdin.
    """
    parse2_py(
        sam_path, chroms_path, output, output_parsed_alignments, output_stats, **kwargs
    )


def parse2_py(
    sam_path, chroms_path, output, output_parsed_alignments, output_stats, **kwargs
):
    ### Set up input stream
    if sam_path:  # open input sam file with pysam
        input_sam = AlignmentFilePairtoolized(
            sam_path, "r", threads=kwargs.get("nproc_in")
        )
    else:  # read from stdin
        input_sam = AlignmentFilePairtoolized("-", "r", threads=kwargs.get("nproc_in"))

    ### Set up output streams
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
    out_alignments_stream = (
        fileio.auto_open(
            output_parsed_alignments,
            mode="w",
            nproc=kwargs.get("nproc_out"),
            command=kwargs.get("cmd_out", None),
        )
        if output_parsed_alignments
        else None
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

    if out_alignments_stream:
        out_alignments_stream.write(
            "readID\tside\tchrom\tpos\tstrand\tmapq\tcigar\tdist_5_lo\tdist_5_hi\tmatched_bp\n"
        )

    # generate empty PairCounter if stats output is requested:
    out_stat = PairCounter() if output_stats else None

    ### Set up output parameters
    add_columns = kwargs.get("add_columns", [])
    add_columns = [col for col in add_columns.split(",") if col]
    for col in add_columns:
        if not (
            (col in pairsam_format.EXTRA_COLUMNS) or (len(col) == 2 and col.isupper())
        ):
            raise Exception("{} is not a valid extra column".format(col))

    columns = pairsam_format.COLUMNS + (
        [c + side for c in add_columns for side in ["1", "2"]]
    )

    if kwargs.get("drop_sam", True):
        columns.pop(columns.index("sam1"))
        columns.pop(columns.index("sam2"))

    if not kwargs.get("add_pair_index", False):
        columns.pop(columns.index("walk_pair_index"))
        columns.pop(columns.index("walk_pair_type"))

    ### Parse header
    samheader = input_sam.header

    if not samheader:
        raise ValueError(
            "The input sam is missing a header! If reading a bam file, please use `samtools view -h` to include the header."
        )

    ### Parse chromosome files present in the input
    sam_chromsizes = headerops.get_chromsizes_from_pysam_header(samheader)
    chromosomes = headerops.get_chrom_order(chroms_path, list(sam_chromsizes.keys()))

    ### Write new header to the pairsam file
    header = headerops.make_standard_pairsheader(
        assembly=kwargs.get("assembly", ""),
        chromsizes=[(chrom, sam_chromsizes[chrom]) for chrom in chromosomes],
        columns=columns,
        shape="whole matrix" if not kwargs["flip"] else "upper triangle",
    )

    header = headerops.insert_samheader_pysam(header, samheader)
    header = headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    outstream.writelines((l + "\n" for l in header))

    ### Parse input and write to the outputs
    streaming_classify(
        input_sam,
        outstream,
        chromosomes,
        out_alignments_stream,
        out_stat,
        parse2=True,
        **kwargs
    )

    # save statistics to a file if it was requested:
    if out_stat:
        out_stat.save(out_stats_stream)

    if outstream != sys.stdout:
        outstream.close()
    if out_alignments_stream:
        out_alignments_stream.close()
    if out_stats_stream:
        out_stats_stream.close()


if __name__ == "__main__":
    parse2()
