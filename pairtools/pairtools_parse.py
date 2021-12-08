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
import pysam

from . import _fileio, _pairsam_format, _parse, _headerops, cli, common_io_options
from .pairtools_stats import PairCounter
from ._parse_pysam import AlignmentFilePairtoolized

UTIL_NAME = "pairtools_parse"

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
    help="output file. "
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4-compressed."
    "By default, the output is printed into stdout. ",
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
    "--max-molecule-size",
    type=int,
    default=750,
    show_default=True,
    help="The maximal size of a Hi-C molecule; used to rescue single ligations"
    "(from molecules with three alignments) and to rescue complex ligations."
    "The default is based on oriented P(s) at short ranges of multiple Hi-C.",
)
@click.option(
    "--drop-readid",
    is_flag=True,
    help="If specified, do not add read ids to the output",
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
    help="If specified, each pair will have junction index in the molecule",
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
    "--output-parsed-alignments",
    type=str,
    default="",
    help="output file for all parsed alignments, including walks."
    " Useful for debugging and rnalysis of walks."
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
@click.option(
    "--report-alignment-end",
    type=click.Choice(["5", "3"]),
    default="5",
    help="specifies whether the 5' or 3' end of the alignment is reported as"
    " the position of the Hi-C read.",
)
@click.option(
    "--max-inter-align-gap",
    type=int,
    default=20,
    show_default=True,
    help="read segments that are not covered by any alignment and"
    ' longer than the specified value are treated as "null" alignments.'
    " These null alignments convert otherwise linear alignments into walks,"
    " and affect how they get reported as a Hi-C pair (see --walks-policy).",
)
@click.option(
    "--walks-policy",
    type=click.Choice(["mask", "5any", "5unique", "3any", "3unique", "all"]),
    default="mask",
    help="the policy for reporting unrescuable walks (reads containing more"
    " than one alignment on one or both sides, that can not be explained by a"
    " single ligation between two mappable DNA fragments)."
    ' "mask" - mask walks (chrom="!", pos=0, strand="-"); '
    ' "5any" - report the 5\'-most alignment on each side;'
    ' "5unique" - report the 5\'-most unique alignment on each side, if present;'
    ' "3any" - report the 3\'-most alignment on each side;'
    ' "3unique" - report the 3\'-most unique alignment on each side, if present;'
    ' "all" - report all available unique alignments on each side.',
    show_default=True,
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
    "--no-flip",
    is_flag=True,
    help="If specified, do not flip pairs in genomic order and instead preserve "
    "the order in which they were sequenced.",
)
@click.option(
    "--pysam-backend",
    is_flag=True,
    help="If specified, parse files with pysam for speedup.",
)
@common_io_options
def parse(
    sam_path,
    chroms_path,
    output,
    assembly,
    min_mapq,
    max_molecule_size,
    drop_readid,
    drop_seq,
    drop_sam,
    add_junction_index,
    add_columns,
    output_parsed_alignments,
    output_stats,
    **kwargs
):
    """Find ligation junctions in .sam, make .pairs.
    SAM_PATH : an input .sam/.bam file with paired-end sequence alignments of
    Hi-C molecules. If the path ends with .bam, the input is decompressed from
    bam with samtools. By default, the input is read from stdin.
    """
    parse_py(
        sam_path,
        chroms_path,
        output,
        assembly,
        min_mapq,
        max_molecule_size,
        drop_readid,
        drop_seq,
        drop_sam,
        add_junction_index,
        add_columns,
        output_parsed_alignments,
        output_stats,
        **kwargs
    )


def parse_py(
    sam_path,
    chroms_path,
    output,
    assembly,
    min_mapq,
    max_molecule_size,
    drop_readid,
    drop_seq,
    drop_sam,
    add_junction_index,
    add_columns,
    output_parsed_alignments,
    output_stats,
    **kwargs
):

    pysam_backend = kwargs.get("pysam_backend", False)

    ### Set up input stream
    if pysam_backend:
        if sam_path:  # open input sam file with pysam
            input_sam = AlignmentFilePairtoolized(sam_path, "r")
        else:  # read from stdin
            input_sam = AlignmentFilePairtoolized("_", "r")
    else:
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

    ### Set up output streams
    outstream = (
        _fileio.auto_open(
            output,
            mode="w",
            nproc=kwargs.get("nproc_out"),
            command=kwargs.get("cmd_out", None),
        )
        if output
        else sys.stdout
    )
    out_alignments_stream = (
        _fileio.auto_open(
            output_parsed_alignments,
            mode="w",
            nproc=kwargs.get("nproc_out"),
            command=kwargs.get("cmd_out", None),
        )
        if output_parsed_alignments
        else None
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

    if out_alignments_stream:
        out_alignments_stream.write(
            "readID\tside\tchrom\tpos\tstrand\tmapq\tcigar\tdist_5_lo\tdist_5_hi\tmatched_bp\n"
        )

    # generate empty PairCounter if stats output is requested:
    out_stat = PairCounter() if output_stats else None

    ### Set up output parameters
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

    ### Parse header
    if pysam_backend:
        samheader = input_sam.header
    else:
        samheader, input_sam = _headerops.get_header(instream, comment_char="@")

    if not samheader:
        raise ValueError(
            "The input sam is missing a header! If reading a bam file, please use `samtools view -h` to include the header."
        )

    ### Parse chromosome files present in the input
    if pysam_backend:
        sam_chromsizes = _headerops.get_chromsizes_from_pysam_header(samheader)
    else:
        sam_chromsizes = _headerops.get_chromsizes_from_sam_header(samheader)
    chromosomes = _headerops.get_chrom_order(chroms_path, list(sam_chromsizes.keys()))

    ### Write new header to the pairsam file
    header = _headerops.make_standard_pairsheader(
        assembly=assembly,
        chromsizes=[(chrom, sam_chromsizes[chrom]) for chrom in chromosomes],
        columns=columns,
        shape="whole matrix" if kwargs["no_flip"] else "upper triangle",
    )

    if pysam_backend:
        header = _headerops.insert_samheader_pysam(header, samheader)
    else:
        header = _headerops.insert_samheader(header, samheader)
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    outstream.writelines((l + "\n" for l in header))

    ### Parse input and write to the outputs
    streaming_classify(
        input_sam,
        outstream,
        chromosomes,
        min_mapq,
        max_molecule_size,
        drop_readid,
        drop_seq,
        drop_sam,
        add_junction_index,
        add_columns,
        out_alignments_stream,
        out_stat,
        **kwargs
    )

    # save statistics to a file if it was requested:
    if out_stat:
        out_stat.save(out_stats_stream)

    if outstream != sys.stdout:
        outstream.close()
    # close optional output streams if needed:
    if out_alignments_stream:
        out_alignments_stream.close()
    if out_stats_stream:
        out_stats_stream.close()


def streaming_classify(
    instream,
    outstream,
    chromosomes,
    min_mapq,
    max_molecule_size,
    drop_readid,
    drop_seq,
    drop_sam,
    add_junction_index,
    add_columns,
    out_alignments_stream,
    out_stat,
    **kwargs
):
    """
    Parse input sam file and write to the outstream(s)
    """

    pysam_backend = kwargs.get("pysam_backend", False)

    ### Store output parameters in usable form:
    chrom_enum = dict(
        zip(
            [_pairsam_format.UNMAPPED_CHROM] + list(chromosomes),
            range(len(chromosomes) + 1),
        )
    )
    sam_tags = [col for col in add_columns if len(col) == 2 and col.isupper()]
    store_seq = "seq" in add_columns

    ### Create temporary variables that will be populated by parsing reads at each iteration over input:
    prev_readID = ""  # Placeholder for the read id
    sams1 = []  # Placeholder for the left alignments
    sams2 = []  # Placeholder for the right alignments
    aligned_segment = ""  # Placeholder for each aligned segment

    ### Compile readID transformation if requested:
    readID_transform = kwargs.get("readid_transform", None)
    if readID_transform is not None:
        readID_transform = compile(readID_transform, "<string>", "eval")

    ### Iterate over the input pysam:
    instream = iter(instream)
    while aligned_segment is not None:
        aligned_segment = next(
            instream, None
        )  # required for proper parsing of the last read

        if pysam_backend:
            readID = aligned_segment.query_name if aligned_segment else None
        else:
            readID = aligned_segment.split("\t", 1)[0] if aligned_segment else None
        if readID_transform is not None and readID is not None:
            readID = eval(readID_transform)

        # Perform parsing and writing when all the segments are parsed from the read:
        if not (aligned_segment) or ((readID != prev_readID) and prev_readID):

            for (
                algn1,
                algn2,
                all_algns1,
                all_algns2,
                junction_index,
            ) in _parse.parse_sams_into_pair(
                sams1,
                sams2,
                min_mapq,
                max_molecule_size,
                kwargs["max_inter_align_gap"],
                kwargs["walks_policy"],
                kwargs["report_alignment_end"] == "3",
                sam_tags,
                store_seq,
                pysam_backend,
            ):

                flip_pair = (not kwargs["no_flip"]) and (
                    not _parse.check_pair_order(algn1, algn2, chrom_enum)
                )

                if flip_pair:
                    algn1, algn2 = algn2, algn1
                    sams1, sams2 = sams2, sams1

                _parse.write_pairsam(
                    algn1,
                    algn2,
                    prev_readID,
                    junction_index,
                    sams1,
                    sams2,
                    outstream,
                    drop_readid,
                    drop_sam,
                    add_junction_index,
                    add_columns,
                    pysam_backend,
                )

                # add a pair to PairCounter if stats output is requested:
                if out_stat:
                    out_stat.add_pair(
                        algn1["chrom"],
                        int(algn1["pos"]),
                        algn1["strand"],
                        algn2["chrom"],
                        int(algn2["pos"]),
                        algn2["strand"],
                        algn1["type"] + algn2["type"],
                    )

                if out_alignments_stream:
                    _parse.write_all_algnments(
                        prev_readID, all_algns1, all_algns2, out_alignments_stream
                    )

            sams1.clear()
            sams2.clear()

        if aligned_segment is not None:
            if pysam_backend:
                _parse.push_pysam(aligned_segment, drop_seq, sams1, sams2)
            else:
                _parse.push_sam(aligned_segment, drop_seq, sams1, sams2)
            prev_readID = readID


if __name__ == "__main__":
    parse()
