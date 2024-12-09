"""
Set of functions used for pairsam parse, migrated from pairtools/pairtools_parse.py

Parse operates with several basic data types:

I. pysam-based:
    1. **sam entry** is a continuous aligned fragment of the read mapped to certain location in the genome.
        Because we read sam entries from .sam/.bam files automatically with modified pysam,
        each sam entry is in fact special AlignedSegmentPairtoolized Cython object
        that has alignment attributes and can be easily accessed from Python.
    
        Sam entries are gathered into reads by `group_alignments_by_side` function.
    
    2. **read** is a collection of sam entries corresponding to a single Hi-C molecule.
        It is represented by three variables:
        readID, sams1 and sams2, which keep left and right sam entries, correspondingly.
        Read is populated from the stream of sam entries on a fly, the process happenning
        in `streaming_classify` function.

II. python-based data types are parsed from pysam-based ones:

    1. **alignment** is a continuous aligned fragment represented as dictionary with relevant fields,
        such as "chrom", "pos5", "pos3", "strand", "type", etc.
        
        `empty_alignment` creates empty alignment,
        `parse_pysam_entry` create new alignments from pysam entries,
        `mask_alignment` clears some fields of the alignment to match the default "unmapped" state.
        
        `flip_alignment`, `flip_orientation` and `flip_ends` are useful functions that help to orient alignments.

    2. **pair** of two alignments is represented by three variables:
        algn1 (left alignment), algn2 (right alignment) and pair_index.
        Pairs are obtained by `parse_read` or `parse2_read`.
        Additionally, these functions also output all alignments for each side.

"""
from . import pairsam_format
from .parse_pysam import get_mismatches_c

def streaming_classify(
    instream, outstream, chromosomes, out_alignments_stream, out_stat, **kwargs
):
    """
    Parse input sam file into individual reads, pairs, walks,
    then write to the outstream(s).

    Additional kwargs:
        min_mapq,
        drop_readid,
        drop_seq,
        drop_sam,
        add_pair_index,
        add_columns, # comma-separated list
        report_alignment_end,
        max_inter_align_gap
    parse:
        max_molecule_size
        walks_policy
    parse2:
        single_end: indicator whether single-end data is provided
        report_position, one of: "outer", "junction", "read", "walk"
        report_orientation, one of: "pair", "junction", "read", "walk"
        dedup_max_mismatch: For intramolecular deduplication
        max_insert_size: maximum insert size when searching for overlapping ends of R1 and R2
        expand: perform combinatorial expansion or not
        max_expansion_depth: maximum expansion depth, works in combination with expand=True

    """

    parse2 = kwargs.get("parse2", False)

    ### Store output parameters in a usable form:
    chrom_enum = dict(
        zip(
            [pairsam_format.UNMAPPED_CHROM] + list(chromosomes),
            range(len(chromosomes) + 1),
        )
    )
    add_columns = kwargs.get("add_columns", "")
    if isinstance(add_columns, str) and len(add_columns) > 0:
        add_columns = add_columns.split(",")
    elif len(add_columns) == 0:
        add_columns = []
    elif not isinstance(add_columns, list):
        raise ValueError(f"Unknown type of add_columns: {type(add_columns)}")

    sam_tags = [col for col in add_columns if len(col) == 2 and col.isupper()]
    store_seq = "seq" in add_columns

    ### Compile readID transformation:
    readID_transform = kwargs.get("readid_transform", None)
    if readID_transform is not None:
        readID_transform = compile(readID_transform, "<string>", "eval")

    ### Iterate over input pysam:
    instream = iter(instream)
    for (readID, (sams1, sams2)) in read_alignment_block(instream, sort=True, group_by_side=True, return_readID=True, readID_transform=readID_transform):

        ### Parse
        if not parse2:  # regular parser:
            pairstream, all_algns1, all_algns2 = parse_read(
                sams1,
                sams2,
                min_mapq=kwargs["min_mapq"],
                max_molecule_size=kwargs["max_molecule_size"],
                max_inter_align_gap=kwargs["max_inter_align_gap"],
                walks_policy=kwargs["walks_policy"],
                sam_tags=sam_tags,
                store_seq=store_seq,
                report_mismatches=True if "mismatches" in add_columns else False,
            )
        else:  # parse2 parser:
            pairstream, all_algns1, all_algns2 = parse2_read(
                sams1,
                sams2,
                min_mapq=kwargs["min_mapq"],
                max_inter_align_gap=kwargs["max_inter_align_gap"],
                max_insert_size=kwargs.get("max_insert_size", 500),
                single_end=kwargs["single_end"],
                report_position=kwargs["report_position"],
                report_orientation=kwargs["report_orientation"],
                sam_tags=sam_tags,
                dedup_max_mismatch=kwargs["dedup_max_mismatch"],
                store_seq=store_seq,
                expand=kwargs["expand"],
                max_expansion_depth=kwargs["max_expansion_depth"],
                report_mismatches=True if "mismatches" in add_columns else False,
            )

        ### Write:
        read_has_alignments = False
        for (algn1, algn2, pair_index) in pairstream:
            read_has_alignments = True

            # Alignment end defaults to 5' if report_alignment_end is unspecified:
            if kwargs.get("report_alignment_end", "5") == "5":
                algn1["pos"] = algn1["pos5"]
                algn2["pos"] = algn2["pos5"]
            else:
                algn1["pos"] = algn1["pos3"]
                algn2["pos"] = algn2["pos3"]

            if kwargs["flip"]:
                flip_pair = not check_pair_order(algn1, algn2, chrom_enum)
                if flip_pair:
                    algn1, algn2 = algn2, algn1
                    sams1, sams2 = sams2, sams1

            write_pairsam(
                algn1,
                algn2,
                readID=readID,
                pair_index=pair_index,
                sams1=sams1,
                sams2=sams2,
                out_file=outstream,
                drop_readid=kwargs["drop_readid"],
                drop_seq=kwargs["drop_seq"],
                drop_sam=kwargs["drop_sam"],
                add_pair_index=kwargs["add_pair_index"],
                add_columns=add_columns,
            )

            # add a pair to PairCounter for stats output:
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

        # write all alignments:
        if out_alignments_stream and read_has_alignments:
            write_all_algnments(
                readID, all_algns1, all_algns2, out_alignments_stream
            )


############################
### Alignment utilities: ###
############################

def empty_alignment():
    return {
        "chrom": pairsam_format.UNMAPPED_CHROM,
        "pos5": pairsam_format.UNMAPPED_POS,
        "pos3": pairsam_format.UNMAPPED_POS,
        "pos": pairsam_format.UNMAPPED_POS,
        "strand": pairsam_format.UNMAPPED_STRAND,
        "dist_to_5": 0,
        "dist_to_3": 0,
        "mapq": 0,
        "is_unique": False,
        "is_mapped": False,
        "is_linear": True,
        "cigar": "*",
        "algn_ref_span": 0,
        "algn_read_span": 0,
        "matched_bp": 0,
        "clip3_ref": 0,
        "clip5_ref": 0,
        "read_len": 0,
        "type": "N",
        "mismatches": "",
    }

def group_alignments_by_side(sams):
    """Group pysam AlignedSegments (sams) into left-read (R1) and right-read (R2) sam entries"""

    sams1 = []
    sams2 = []
    for sam_entry in sams:
        flag = sam_entry.flag
        if (flag & 0x40) != 0:
            sams1.append(sam_entry)  # left read, or first read in a pair
        else:
            sams2.append(sam_entry)  # right read, or mate pair
    return sams1, sams2


def read_alignment_block(instream, sort=True, group_by_side=True, return_readID=True, readID_transform=None):
    sams = []

    prev_readID = None
    while True:
        sam_entry = next(instream, None)
        readID = sam_entry.query_name if sam_entry else None
        if readID_transform is not None and readID is not None: 
            readID = eval(readID_transform)

        # Read is fully populated, then parse and write:
        if not (sam_entry) or ((readID != prev_readID) and prev_readID):
            if sort:
                sams = sorted(sams, key=lambda a: (a.is_read2, a.query_alignment_start))
            out = sams if not group_by_side else group_alignments_by_side(sams)
            out = out if not return_readID else (prev_readID, out)
            yield out
            
            sams.clear()
            
        if sam_entry is None:
            break
        else:
            sams.append(sam_entry)
            prev_readID = readID

def parse_pysam_entry(
    sam,
    min_mapq,
    sam_tags=None,
    store_seq=False,
    report_3_alignment_end=False,
    report_mismatches=False,
):
    """Parse alignments from pysam AlignedSegment entry
    :param sam: input pysam AlignedSegment entry
    :param min_mapq: minimal MAPQ to consider as a proper alignment
    :param sam_tags: list of sam tags to store
    :param store_seq: if True, the sequence will be parsed and stored in the output
    :param report_3_alignment_end: if True, 3'-end of alignment will be
                                reported as position (will be deprecated)
    :param report_mismatches: if True, mismatches will be parsed from MD field
    :return: parsed aligned entry (dictionary)
    """

    flag = sam.flag
    is_mapped = (flag & 0x04) == 0
    mapq = sam.mapq
    is_unique = sam.is_unique(min_mapq)
    is_linear = sam.is_linear
    cigar = sam.cigar_dict
    if is_mapped:
        if (flag & 0x10) == 0:
            strand = "+"
            dist_to_5 = cigar["clip5_ref"]
            dist_to_3 = cigar["clip3_ref"]
        else:
            strand = "-"
            dist_to_5 = cigar["clip3_ref"]
            dist_to_3 = cigar["clip5_ref"]

        if is_unique:
            chrom = sam.reference_name
            if strand == "+":
                # Note that pysam output is zero-based, thus add +1:
                pos5 = sam.reference_start + 1
                pos3 = sam.reference_start + cigar["algn_ref_span"]
            else:
                pos5 = sam.reference_start + cigar["algn_ref_span"]
                # Note that pysam output is zero-based, thus add +1:
                pos3 = sam.reference_start + 1

            # Get number of matches:
            if not sam.has_tag("MD") or not report_mismatches:
                mismatches = ""
            else:
                seq = sam.query_sequence.upper()
                quals = sam.query_qualities
                aligned_pairs = sam.get_aligned_pairs(with_seq=True, matches_only=True)
                mismatches = get_mismatches_c(seq, quals, aligned_pairs)
                mismatches = ",".join(
                    [
                        f"{original}:{mutated}:{phred}:{ref}:{read}"
                        for original, mutated, phred, ref, read in mismatches
                    ]
                )
                # n_matches = len(aligned_pairs)

        else:
            chrom = pairsam_format.UNMAPPED_CHROM
            strand = pairsam_format.UNMAPPED_STRAND
            pos5 = pairsam_format.UNMAPPED_POS
            pos3 = pairsam_format.UNMAPPED_POS
            mismatches = ""
    else:
        chrom = pairsam_format.UNMAPPED_CHROM
        strand = pairsam_format.UNMAPPED_STRAND
        pos5 = pairsam_format.UNMAPPED_POS
        pos3 = pairsam_format.UNMAPPED_POS

        dist_to_5 = 0
        dist_to_3 = 0
        mismatches = ""

    algn = {
        "chrom": chrom,
        "pos5": pos5,
        "pos3": pos3,
        "strand": strand,
        "mapq": mapq,
        "is_mapped": is_mapped,
        "is_unique": is_unique,
        "is_linear": is_linear,
        "dist_to_5": dist_to_5,
        "dist_to_3": dist_to_3,
        "type": ("N" if not is_mapped else ("M" if not is_unique else "U")),
        "mismatches": mismatches,
    }

    algn.update(cigar)

    algn["pos"] = algn["pos3"] if report_3_alignment_end else algn["pos5"]

    ### Add tags to the alignment:
    if sam_tags:
        tags = sam.tags
        for tag in sam_tags:
            algn[tag] = ""

        for col, value in tags:
            for tag in sam_tags:
                if col == tag:
                    algn[tag] = value
                    continue

    if store_seq:
        algn["seq"] = sam.seq

    return algn


def mask_alignment(algn):
    """
    Reset the coordinates of an alignment.
    """
    algn["chrom"] = pairsam_format.UNMAPPED_CHROM
    algn["pos5"] = pairsam_format.UNMAPPED_POS
    algn["pos3"] = pairsam_format.UNMAPPED_POS
    algn["pos"] = pairsam_format.UNMAPPED_POS
    algn["strand"] = pairsam_format.UNMAPPED_STRAND

    return algn


def flip_alignment(hic_algn):
    """
    Flip a single alignment as if it was sequenced from the opposite end
    :param hic_algn: Alignment to be modified
    :return:
    """
    hic_algn = dict(hic_algn)  # overwrite the variable with the copy of dictionary
    hic_algn["pos5"], hic_algn["pos3"] = hic_algn["pos3"], hic_algn["pos5"]
    hic_algn["strand"] = "+" if (hic_algn["strand"] == "-") else "-"
    return hic_algn


def flip_orientation(hic_algn):
    """
    Flip orientation of a single alignment
    :param hic_algn: Alignment to be modified
    :return:
    """
    hic_algn = dict(hic_algn)  # overwrite the variable with the copy of dictionary
    hic_algn["strand"] = "+" if (hic_algn["strand"] == "-") else "-"
    return hic_algn


def flip_position(hic_algn):
    """
    Flip ends of a single alignment
    :param hic_algn: Alignment to be modified
    :return:
    """
    hic_algn = dict(hic_algn)  # overwrite the variable with the copy of dictionary
    hic_algn["pos5"], hic_algn["pos3"] = hic_algn["pos3"], hic_algn["pos5"]
    return hic_algn



def _convert_gaps_into_alignments(sorted_algns, max_inter_align_gap):
    """
    Inplace conversion of gaps longer than max_inter_align_gap into alignments
    """
    if (len(sorted_algns) == 1) and (not sorted_algns[0]["is_mapped"]):
        return

    last_5_pos = 0
    for i in range(len(sorted_algns)):
        algn = sorted_algns[i]
        if algn["dist_to_5"] - last_5_pos > max_inter_align_gap:
            new_algn = empty_alignment()
            new_algn["dist_to_5"] = last_5_pos
            new_algn["algn_read_span"] = algn["dist_to_5"] - last_5_pos
            new_algn["read_len"] = algn["read_len"]
            new_algn["dist_to_3"] = new_algn["read_len"] - algn["dist_to_5"]

            last_5_pos = algn["dist_to_5"] + algn["algn_read_span"]

            sorted_algns.insert(i, new_algn)
            i += 2
        else:
            last_5_pos = max(last_5_pos, algn["dist_to_5"] + algn["algn_read_span"])
            i += 1


def normalize_alignment_list(algns, side, sort_by="dist_to_5", max_inter_align_gap=None):
    """
    Normalize the alignment list: insert empty alignments in gaps between alignments, 
    sort by distance to the 5' end, add read side, alignment index.

    Args:
        algns (list): The list of alignments.
        side (str): The side of the alignment.
        sort_by (str, optional): The key to sort the alignments by. Defaults to "dist_to_5".
        max_inter_align_gap (int, optional): The maximum allowed gap between alignments. Defaults to None.

    Returns:
        list: The normalized alignment list.

    """
    if len(algns) == 0:
        algns = [empty_alignment()]

    if sort_by:
        algns = sorted(algns, key=lambda algn: algn[sort_by])

    if max_inter_align_gap is not None:
        _convert_gaps_into_alignments(algns, max_inter_align_gap)

    for i, algn in enumerate(algns):
        algn["read_side"] = side
        algn["algn_idx"] = i
        algn["same_side_algn_count"] = len(algns)

    return algns
            


####################
### Parsing utilities:
####################


def parse_read(
    sams1,
    sams2,
    min_mapq,
    max_molecule_size,
    max_inter_align_gap,
    walks_policy,
    sam_tags,
    store_seq,
    report_mismatches=False,
):
    """
    Parse sam entries corresponding to a single read (or Hi-C molecule)
    into pairs of alignments.

    Returns
    -------
    stream: iterator
        Each element is a triplet: (algn1, aldn2, pair_index)
        algn1, algn2: dict
            Two alignments selected for reporting as a Hi-C pair.
        pair_index
            pair index of a pair in the molecule.
    algns1, algns2: lists
        All alignments, sorted according to their order in on a read.
    """

    # Check if there is at least one sam entry per side:
    if walks_policy == "all":
        is_empty = (len(sams1) == 0 and len(sams2) < 2) or (
            len(sams2) == 0 and len(sams1) < 2
        )
    else:
        is_empty = (len(sams1) == 0) or (len(sams2) == 0)
    if is_empty:
        algns1 = [empty_alignment()]
        algns2 = [empty_alignment()]
        algns1[0]["type"] = "X"
        algns2[0]["type"] = "X"
        pair_index = (1, "R1-2")
        return iter([(algns1[0], algns2[0], pair_index)]), algns1, algns2

    # Generate a sorted, gap-filled list of all alignments
    algns1 = [
        parse_pysam_entry(
            sam, min_mapq, sam_tags, store_seq, report_mismatches=report_mismatches
        )
        for sam in sams1
    ]
    algns2 = [
        parse_pysam_entry(
            sam, min_mapq, sam_tags, store_seq, report_mismatches=report_mismatches
        )
        for sam in sams2
    ]

    algns1 = normalize_alignment_list(algns1, 1, sort_by="dist_to_5", max_inter_align_gap=max_inter_align_gap)
    algns2 = normalize_alignment_list(algns2, 2, sort_by="dist_to_5", max_inter_align_gap=max_inter_align_gap)


    # By default, assume each molecule is a single pair with single unconfirmed pair:
    hic_algn1 = algns1[0]
    hic_algn2 = algns2[0]
    pair_index = (1, "R1-2")

    # Define the type of alignment on each side:
    is_chimeric_1 = len(algns1) > 1
    is_chimeric_2 = len(algns2) > 1

    # Parse chimeras
    if is_chimeric_1 or is_chimeric_2:

        # Report all the linear alignments in a read pair
        if walks_policy == "all":
            # Report linear alignments after deduplication of complex walks with default settings:
            return (
                parse_complex_walk(
                    algns1,
                    algns2,
                    max_molecule_size,
                    report_position="outer",
                    report_orientation="pair",
                ),
                algns1,
                algns2,
            )

        elif walks_policy in ["mask", "5any", "5unique", "3any", "3unique"]:
            # Report only two alignments for a read pair
            rescued_linear_side = rescue_walk(algns1, algns2, max_molecule_size)

            # Walk was rescued as a simple walk:
            if rescued_linear_side is not None:
                pair_index = (1, "R1" if rescued_linear_side == 1 else "R2")
            # Walk is unrescuable:
            else:
                if walks_policy == "mask":
                    hic_algn1 = mask_alignment(dict(hic_algn1))
                    hic_algn2 = mask_alignment(dict(hic_algn2))
                    hic_algn1["type"] = "W"
                    hic_algn2["type"] = "W"

                elif walks_policy == "5any":
                    hic_algn1 = algns1[0]
                    hic_algn2 = algns2[0]

                elif walks_policy == "5unique":
                    hic_algn1 = algns1[0]
                    for algn in algns1:
                        if algn["is_mapped"] and algn["is_unique"]:
                            hic_algn1 = algn
                            break

                    hic_algn2 = algns2[0]
                    for algn in algns2:
                        if algn["is_mapped"] and algn["is_unique"]:
                            hic_algn2 = algn
                            break

                elif walks_policy == "3any":
                    hic_algn1 = algns1[-1]
                    hic_algn2 = algns2[-1]

                elif walks_policy == "3unique":
                    hic_algn1 = algns1[-1]
                    for algn in algns1[::-1]:
                        if algn["is_mapped"] and algn["is_unique"]:
                            hic_algn1 = algn
                            break

                    hic_algn2 = algns2[-1]
                    for algn in algns2[::-1]:
                        if algn["is_mapped"] and algn["is_unique"]:
                            hic_algn2 = algn
                            break

                # DEPRECATED: lower-case reported walks on the chimeric side
                if walks_policy != "mask":
                    if is_chimeric_1:
                        hic_algn1 = dict(hic_algn1)
                        # hic_algn1["type"] = hic_algn1["type"].lower()
                    if is_chimeric_2:
                        hic_algn2 = dict(hic_algn2)
                        # hic_algn2["type"] = hic_algn2["type"].lower()

        else:
            raise ValueError(f"Walks policy {walks_policy} is not supported.")

    return iter([(hic_algn1, hic_algn2, pair_index)]), algns1, algns2


def parse2_read(
    sams1,
    sams2,
    min_mapq,
    max_inter_align_gap,
    max_insert_size,
    single_end,
    report_position="outer",
    report_orientation="pair",
    sam_tags=[],
    dedup_max_mismatch=3,
    store_seq=False,
    report_mismatches=False,
    expand=False,
    max_expansion_depth=None,
):
    """
    Parse sam entries corresponding to a Hi-C molecule into alignments in parse2 mode
    for a Hi-C pair.
    Returns
    -------
    stream: iterator
        Each element is a triplet: (algn1, aldn2, pair_index)
        algn1, algn2: dict
            Two alignments selected for reporting as a Hi-C pair.
        pair_index
            pair index of a pair in the molecule, a tuple: (1, "R1-2")
    algns1, algns2: lists
        All alignments, sorted according to their order in on a read.
    """

    # Single-end mode:
    if single_end:
        # Generate a sorted, gap-filled list of all alignments
        algns1 = [
            parse_pysam_entry(
                sam, min_mapq, sam_tags, store_seq, report_mismatches=report_mismatches
            )
            for sam in sams2  # note sams2, that's how these reads are typically parsed
        ]
        
        algns1 = normalize_alignment_list(algns1, 1, sort_by="dist_to_5", max_inter_align_gap=max_inter_align_gap)
        algns2 = []  # Empty alignment dummy

        if len(algns1) > 1:
            # Look for ligation pair, and report linear alignments after deduplication of complex walks:
            # (Note that coordinate system for single-end reads does not change the behavior)
            output = parse_complex_walk(
                algns1,
                algns2,
                max_insert_size,
                report_position,
                report_orientation,
                dedup_max_mismatch,
                expand,
                max_expansion_depth,
            )
            output = [x for x in output if x[-1][-1] != "R1-2"]
            return (output, algns1, algns2)
        elif len(algns1) == 1:
            # If no additional information, we assume each molecule is a single ligation with single unconfirmed pair:
            algn2 = empty_alignment()
            pair_index = (1, "R1")
            return iter([(algns1[0], algn2, pair_index)]), algns1, algns2
        else:
            # If no additional information, we assume each molecule is a single ligation with single unconfirmed pair:
            algn1 = empty_alignment()
            algn2 = empty_alignment()
            pair_index = (1, "R1")
            return iter([(algn1, algn2, pair_index)]), algns1, algns2

    # Paired-end mode:
    else:
        # Check if there is at least one SAM entry per side:
        is_empty = (len(sams1) == 0 and len(sams2) < 2) or (
            len(sams2) == 0 and len(sams1) < 2
        )
        if is_empty:
            algns1 = [empty_alignment()]
            algns2 = [empty_alignment()]
            algns1[0]["type"] = "X"
            algns2[0]["type"] = "X"
            pair_index = (1, "R1-2")
            return iter([(algns1[0], algns2[0], pair_index)]), algns1, algns2

        # Generate a sorted, gap-filled list of all alignments
        algns1 = [
            parse_pysam_entry(
                sam, min_mapq, sam_tags, store_seq, report_mismatches=report_mismatches
            )
            for sam in sams1
        ]
        algns2 = [
            parse_pysam_entry(
                sam, min_mapq, sam_tags, store_seq, report_mismatches=report_mismatches
            )
            for sam in sams2
        ]

        algns1 = normalize_alignment_list(algns1, 1, sort_by="dist_to_5", max_inter_align_gap=max_inter_align_gap)
        algns2 = normalize_alignment_list(algns2, 2, sort_by="dist_to_5", max_inter_align_gap=max_inter_align_gap)

        is_chimeric_1 = len(algns1) > 1
        is_chimeric_2 = len(algns2) > 1

        if is_chimeric_1 or is_chimeric_2:
            # If at least one side is chimera, we must look for ligation pair, and
            # report linear alignments after deduplication of complex walks:
            return (
                parse_complex_walk(
                    algns1,
                    algns2,
                    max_insert_size,
                    report_position,
                    report_orientation,
                    dedup_max_mismatch,
                    expand,
                    max_expansion_depth,
                ),
                algns1,
                algns2,
            )
        else:
            # If no additional information, we assume each molecule is a single ligation with single unconfirmed pair:
            algn2 = algns2[0]
            if report_orientation == "walk":
                algn2 = flip_orientation(algn2)
            if report_position == "walk":
                algn2 = flip_position(algn2)
            pair_index = (1, "R1-2")
            return iter([(algns1[0], algn2, pair_index)]), algns1, algns2


####################
### Walks utilities:
####################


def rescue_walk(algns1, algns2, max_molecule_size):
    """
    Rescue a single ligation that appears as a walk.
    Checks if a molecule with three alignments could be formed via a single
    ligation between two fragments, where one fragment was so long that it
    got sequenced on both sides.
    Uses three criteria:
    1) the 3'-end alignment on one side maps to the same chromosome as the
    alignment fully covering the other side (i.e. the linear alignment)
    2) the two alignments point towards each other on the chromosome
    3) the distance between the outer ends of the two alignments is below
    the specified threshold.
    Alternatively, a single ligation get rescued when the 3' sub-alignment
    maps to multiple locations or no locations at all.

    In the case of a successful rescue, tags the 3' sub-alignment with
    type='X' and the linear alignment on the other side with type='R'.
    Returns
    -------
    linear_side : int
        If the case of a successful rescue, returns the index of the side
        with a linear alignment.
    """

    # If both sides have one alignment or none, no need to rescue!
    n_algns1 = len(algns1)
    n_algns2 = len(algns2)

    if (n_algns1 <= 1) and (n_algns2 <= 1):
        return None

    # Can rescue only pairs with one chimeric alignment with two parts.
    if not (
        ((n_algns1 == 1) and (n_algns2 == 2)) or ((n_algns1 == 2) and (n_algns2 == 1))
    ):
        return None

    first_read_is_chimeric = n_algns1 > 1
    chim5_algn = algns1[0] if first_read_is_chimeric else algns2[0]
    chim3_algn = algns1[1] if first_read_is_chimeric else algns2[1]
    linear_algn = algns2[0] if first_read_is_chimeric else algns1[0]

    # the linear alignment must be uniquely mapped
    if not (linear_algn["is_mapped"] and linear_algn["is_unique"]):
        return None

    can_rescue = True
    # we automatically rescue chimeric alignments with null and non-unique
    # alignments at the 3' side
    if chim3_algn["is_mapped"] and chim5_algn["is_unique"]:
        # 1) in rescued walks, the 3' alignment of the chimeric alignment must be on
        # the same chromosome as the linear alignment on the opposite side of the
        # molecule
        can_rescue &= chim3_algn["chrom"] == linear_algn["chrom"]

        # 2) in rescued walks, the 3' supplemental alignment of the chimeric
        # alignment and the linear alignment on the opposite side must point
        # towards each other
        can_rescue &= chim3_algn["strand"] != linear_algn["strand"]
        if linear_algn["strand"] == "+":
            can_rescue &= linear_algn["pos5"] < chim3_algn["pos5"]
        else:
            can_rescue &= linear_algn["pos5"] > chim3_algn["pos5"]

        # 3) in single ligations appearing as walks, we can infer the size of
        # the molecule and this size must be smaller than the maximal size of
        # Hi-C molecules after the size selection step of the Hi-C protocol
        if linear_algn["strand"] == "+":
            molecule_size = (
                chim3_algn["pos5"]
                - linear_algn["pos5"]
                + chim3_algn["dist_to_5"]
                + linear_algn["dist_to_5"]
            )
        else:
            molecule_size = (
                linear_algn["pos5"]
                - chim3_algn["pos5"]
                + chim3_algn["dist_to_5"]
                + linear_algn["dist_to_5"]
            )

        can_rescue &= molecule_size <= max_molecule_size

    if can_rescue:
        # changing the type of the 3' alignment on side 1, does not show up in the output:
        if first_read_is_chimeric:

            algns1[1]["type"] = "X"
            algns2[0]["type"] = "R"
            return 1
        # changing the type of the 3' alignment on side 2, does not show up in the output:
        else:
            algns1[0]["type"] = "R"
            algns2[1]["type"] = "X"
            return 2
    else:
        return None


def parse_complex_walk(
    algns1,
    algns2,
    max_insert_size,
    report_position,
    report_orientation,
    dedup_max_mismatch=3,
    expand=False,
    max_expansion_depth=None,
):
    """
    Parse a set of ligations that appear as a complex walk.
    This procedure is equivalent to intramolecular deduplication that preserved pair order in a walk.

    :param algns1: List of sequential lefts alignments
    :param algns2: List of sequential right alignments
    :param max_insert_size: maximum insert size when searching for overlapping ends of R1 and R2
    :param report_position: one of "outer", "junction", "read", "walk"; sets pos5 and pos3
    :param report_orientation: one of "pair", "junction", "read", "walk"; sets strand
    :param dedup_max_mismatch: allowed mismatch between intramolecular alignments to detect readthrough duplicates
    :param expand: perform combinatorial expansion of pairs or not
    :param max_expansion_depth: maximum depth (number of segments separating pair). All by default.

    :return: iterator with parsed pairs

    **Intramolecular deduplication**

    Forward read (left):                       right read (right):
    5'------------------------->3'     3'<--------------------------5'
            algns1                              algns2
    <5---3><5---3><5---3><5---3>        <3---5><3---5><3---5><3---5>
        l0     l1    l2     l3              r3     r2     r1    r0

    Alignment - bwa mem reported hit or alignment after gaps conversion.
    Left and right alignments (algns1: [l0, l1, l2, l3], algns2: [r0, r1, r2, r3])
    - alignments on left and right reads reported from 5' to 3' orientation.

    Intramolecular deduplication consists of two steps:
    I. iterative search of overlapping alignment pairs (aka overlap),
    II. if no overlaps or search not possible (less than 2 alignments on either sides),
    search for overlap of end alignments (aka partial overlap).
    III. report pairs before the overlap, deduplicated pairs of overlap and pairs after that.

    Iterative search of overlap is in fact scanning of the right read pairs for the hit
    with the 3'-most pair of the left read:
        1. Initialize.
            Start from 3' of left and right reads. Set `current_left_pair` and `current_right_pair` pointers
        2. Initial compare.
            Compare pairs l2-l3 and r3-r2 by `pairs_overlap`.
                If successful, we found the overlap, go to reporting.
                If unsuccessful, continue search.
        3. Increment.
            Shift `current_right_pair` pointer by one (e.g., take the pair r2-r1).
        4. Check.
            Check that this pair can form a potential overlap with left alignments:
            the number of pairs downstream from l2-l3 on left read should not be less than
            the number of pairs upstream from r2-r1 on right read.
                If overlap cannot be formed, no other overlap in this complex walk is possible, safely exit.
                If the potential overlap can be formed, continue comparison.
        5. Compare.
            Compare the current pair of pairs on left and right reads.
                If comparison fails, go to step 3.
                If comparison is successful, go to 6.
        6. Verify.
            Check that downstream pairs on the left read overlap with the upstream pairs on the right read.
                If yes, exit.
                If not, we do not have an overlap, go to step 3.
    """

    AVAILABLE_REPORT_POSITION = ["outer", "junction", "read", "walk"]
    assert report_position in AVAILABLE_REPORT_POSITION, (
        f"Cannot report position {report_position}, as it is not implemented"
        f'Available choices are: {", ".join(AVAILABLE_REPORT_POSITION)}'
    )

    AVAILABLE_REPORT_ORIENTATION = ["pair", "junction", "read", "walk"]
    assert report_orientation in AVAILABLE_REPORT_ORIENTATION, (
        f"Cannot report orientation {report_orientation}, as it is not implemented"
        f'Available choices are: {", ".join(AVAILABLE_REPORT_ORIENTATION)}'
    )

    output_pairs = []

    # Initialize (step 1).
    n_algns1 = len(algns1)
    n_algns2 = len(algns2)
    current_left_pair = current_right_pair = 1
    remaining_left_pairs = (
        n_algns1 - 1
    )  # Number of possible pairs remaining on left read
    remaining_right_pairs = (
        n_algns2 - 1
    )  # Number of possible pairs remaining on right read
    checked_right_pairs = (
        0  # Number of checked pairs on right read (from the end of read)
    )
    is_overlap = False

    # I. Iterative search of overlap, at least two alignments on each side:
    if (n_algns1 >= 2) and (n_algns2 >= 2):
        # Iteration includes check (step 4):
        while (remaining_left_pairs > checked_right_pairs) and (
            remaining_right_pairs > 0
        ):
            pair1 = (algns1[-current_left_pair - 1], algns1[-current_left_pair])
            pair2 = (algns2[-current_right_pair - 1], algns2[-current_right_pair])
            # Compare (initial or not, step 2 or 5):
            is_overlap = pairs_overlap(
                pair1, pair2, dedup_max_mismatch=dedup_max_mismatch
            )
            if is_overlap:
                last_idx_left_temp = current_left_pair
                last_idx_right_temp = current_right_pair
                checked_right_temp = checked_right_pairs
                # Verify (step 6):
                while is_overlap and (checked_right_temp > 0):
                    last_idx_left_temp += 1
                    last_idx_right_temp -= 1
                    pair1 = (
                        algns1[-last_idx_left_temp - 1],
                        algns1[-last_idx_left_temp],
                    )
                    pair2 = (
                        algns2[-last_idx_right_temp - 1],
                        algns2[-last_idx_right_temp],
                    )
                    is_overlap &= pairs_overlap(
                        pair1, pair2, dedup_max_mismatch=dedup_max_mismatch
                    )
                    checked_right_temp -= 1
                if is_overlap:  # exit
                    current_right_pair += 1
                    break

            # Increment pointers (step 3)
            current_right_pair += 1
            checked_right_pairs += 1
            remaining_right_pairs -= 1

        # No overlap found, roll the current_idx_right back to the initial value:
        if not is_overlap:
            current_right_pair = 1

    if (n_algns2 == 0):
        last_reported_alignment_left = 1
        last_reported_alignment_right = 0
    else:
        # II. Search of partial overlap if there are less than 2 alignments at either sides, or no overlaps found
        if (current_right_pair == 1):
            last_reported_alignment_left = last_reported_alignment_right = 1
            if partial_overlap(
                algns1[-1],
                algns2[-1],
                max_insert_size=max_insert_size,
                dedup_max_mismatch=dedup_max_mismatch,
            ):
                if (
                    n_algns1 >= 2
                ):  # single alignment on right read and multiple alignments on left
                    pair_index = (len(algns1) - 1, "R1")
                    output_pairs.append(
                        format_pair(
                            algns1[-2],
                            algns1[-1],
                            pair_index=pair_index,
                            algn2_pos3=algns2[-1]["pos5"],
                            report_position=report_position,
                            report_orientation=report_orientation,
                        )
                    )
                    last_reported_alignment_left = 2  # set the pointer for reporting

                if (
                    n_algns2 >= 2
                ):  # single alignment on left read and multiple alignments on right
                    pair_index = (len(algns1), "R2")
                    output_pairs.append(
                        format_pair(
                            algns2[-1],
                            algns2[-2],
                            pair_index=pair_index,
                            algn1_pos3=algns1[-1]["pos5"],
                            report_position=report_position,
                            report_orientation=report_orientation,
                        )
                    )
                    last_reported_alignment_right = 2  # set the pointer for reporting

                # Note that if n_algns1==n_algns2==1 and alignments overlap, then we don't need to check,
                # it's a non-ligated DNA fragment that we don't report.

            else:  # end alignments do not overlap, report regular pair:
                pair_index = (len(algns1), "R1-2")
                output_pairs.append(
                    format_pair(
                        algns1[-1],
                        algns2[-1],
                        pair_index=pair_index,
                        report_position=report_position,
                        report_orientation=report_orientation,
                    )
                )

        else:  # there was an overlap, set some pointers:
            last_reported_alignment_left = (
                last_reported_alignment_right
            ) = current_right_pair

    # III. Report all remaining alignments.
    # Report all unique alignments on left read (sequential):
    for i in range(0, n_algns1 - last_reported_alignment_left):
        pair_index = (i + 1, "R1")
        output_pairs.append(
            format_pair(
                algns1[i],
                algns1[i + 1],
                pair_index=pair_index,
                report_position=report_position,
                report_orientation=report_orientation,
            )
        )

    # Report the pairs where both left alignments overlap right:
    for i_overlapping in range(current_right_pair - 1):
        idx_left = n_algns1 - current_right_pair + i_overlapping
        idx_right = n_algns2 - 1 - i_overlapping
        pair_index = (idx_left + 1, "R1&2")
        output_pairs.append(
            format_pair(
                algns1[idx_left],
                algns1[idx_left + 1],
                pair_index=pair_index,
                algn2_pos3=algns2[idx_right - 1]["pos5"],
                report_position=report_position,
                report_orientation=report_orientation,
            )
        )

    # Report all the sequential chimeric pairs in the right read, but not the overlap:
    reporting_order = range(
        0, min(current_right_pair, n_algns2 - last_reported_alignment_right)
    )
    for i in reporting_order:
        # Determine the pair index depending on what is the overlap:
        shift = -1 if current_right_pair > 1 else 0
        pair_index = (
            (
                n_algns1
                + min(current_right_pair, n_algns2 - last_reported_alignment_right)
                - i
                + shift
            ),
            "R2",
        )
        output_pairs.append(
            format_pair(
                algns2[i + 1],
                algns2[i],
                pair_index=pair_index,
                report_position=report_position,
                report_orientation=report_orientation,
            )
        )

    # Sort the pairs according to the pair index:
    output_pairs.sort(key=lambda x: int(x[-1][0]))
    if expand:
        output_pairs = expand_pairs(output_pairs, max_expansion_depth)
    return iter(output_pairs)


### Additional functions for pairs ###
def expand_pairs(pairs_list, max_expansion_depth=None):
    """
    Perform combinatorial expansion of the pairs.

    Parameters
    ----------
    pairs_list: List of formatted pairs (triplets: algn1, algn2, pair_index).
    max_expansion_depth: maximum depth of expansion; all by default (None),
                    number will enforce only pairs from the same strand.
    Returns
    -------
    list of expanded pairs

    """
    for algn1, _algn1, pair_index1 in pairs_list:
        for _algn2, algn2, pair_index2 in pairs_list:
            if pair_index1 > pair_index2:
                continue
            elif pair_index1 == pair_index2:
                # output regular pair with no change
                yield algn1, _algn1, pair_index1
            else:
                pair_order1, pair_type1 = pair_index1
                pair_order2, pair_type2 = pair_index2
                separated_by = pair_order2 - pair_order1
                if (
                    pair_type1 == "R1-2"
                    or pair_type2 == "R1-2"
                    or (pair_type1 == "R1" and pair_type2 == "R2")
                ):
                    pair_type = "R1-2"
                elif pair_type1 == pair_type2:
                    pair_type = pair_type1
                elif pair_type1 == "R1&2":
                    pair_type = pair_type2
                elif pair_type2 == "R1&2":
                    pair_type = pair_type1
                else:
                    raise ValueError(
                        f"Unexpected error, pair types: {pair_type1}, {pair_type2}"
                    )
                same_read = pair_type != "R1-2"

                if (max_expansion_depth is None) or (
                    (separated_by <= max_expansion_depth) and same_read
                ):
                    pair_type = f"E{separated_by}_{pair_type}"
                    yield algn1, algn2, (pair_order1, pair_type)


### Additional functions for complex walks rescue ###
def partial_overlap(algn1, algn2, max_insert_size=500, dedup_max_mismatch=5):
    """
    Two ends of alignments overlap if:
     1) they are from the same chromosome,
     2) map in the opposite directions,
     3) the distance between the outer ends of the two alignments is below the specified max_insert_size,
     4) the distance between the outer ends of the two alignments is above the maximum alignment size.
    (4) guarantees that the alignments point towards each other on the chromosomes.

    Allowed mismatch between intramolecular alignments to detect readthrough duplicates.

    Return: 1 if the alignments overlap or both have troubles with unique mapping,
            0 if they do not overlap or if we do not have enough information
            (e.g. only one of the alignments have troubles with being mapped).
    """

    # Alignments with no match or with multiple matches are counted as overlaps
    if not (algn1["is_mapped"] and algn1["is_unique"]):
        if not (algn2["is_mapped"] and algn2["is_unique"]):
            return 1

    # We assume that successful alignment cannot be an overlap with unmapped or multi-mapped region
    if not (algn1["is_mapped"] and algn1["is_unique"]):
        return 0
    if not (algn2["is_mapped"] and algn2["is_unique"]):
        return 0

    # Both alignments are mapped and unique
    do_overlap = True

    do_overlap &= algn1["chrom"] == algn2["chrom"]
    do_overlap &= algn1["strand"] != algn2["strand"]

    if algn1["strand"] == "+":
        min_algn_size = max(
            algn1["pos3"] - algn1["pos5"], algn2["pos5"] - algn2["pos3"]
        )
        distance_outer_ends = algn2["pos5"] - algn1["pos5"]
    else:
        min_algn_size = max(
            algn1["pos5"] - algn1["pos3"], algn2["pos3"] - algn2["pos5"]
        )
        distance_outer_ends = algn1["pos5"] - algn2["pos5"]

    do_overlap &= distance_outer_ends <= max_insert_size + dedup_max_mismatch
    do_overlap &= distance_outer_ends >= min_algn_size - dedup_max_mismatch

    if do_overlap:
        return 1
    return 0


def pairs_overlap(algns1, algns2, dedup_max_mismatch=3):
    """
    We assume algns1 originate from left read, and algns2 originate from right read:
    left read:                             right read:
    ---------------------------->     <----------------------------
                algns1                             algns2
    5------------3_5------------3     3------------5_3------------5'
    left_5'-algn    left_3'-algn      right_3'-algn   right_5'-algn

    Two pairs of alignments overlap if:
    1) chromosomes/mapping/strand of left_5'-algn and right_3'-algn are the same,
    2) chromosomes/mapping/strand of left_3'-algn and right_5'-algn are the same,
    3) pos3 of left_5'-algn is close to pos5 of right_3'-algn (with dedup_max_mismatch), and
    4) pos5 of left_3'-algn is close to pos3 of right_5'-algn.

    Return: 1 of the pairs of alignments overlap, 0 otherwise.
    """
    left5_algn = algns1[0]
    left3_algn = algns1[1]
    right5_algn = algns2[0]
    right3_algn = algns2[1]

    # We assume that successful alignment cannot be an overlap with unmapped or multi-mapped region:
    mapped_left5_algn = left5_algn["is_mapped"] and left5_algn["is_unique"]
    mapped_left3_algn = left3_algn["is_mapped"] and left3_algn["is_unique"]
    mapped_right5_algn = right5_algn["is_mapped"] and right5_algn["is_unique"]
    mapped_right3_algn = right3_algn["is_mapped"] and right3_algn["is_unique"]

    if not mapped_left5_algn and not mapped_right3_algn:
        left_overlap = True
    elif not mapped_left5_algn and mapped_right3_algn:
        left_overlap = False
    elif mapped_left5_algn and not mapped_right3_algn:
        left_overlap = False
    else:
        left_overlap = True
        left_overlap &= left5_algn["chrom"] == right3_algn["chrom"]
        left_overlap &= left5_algn["strand"] != right3_algn["strand"]

    if not mapped_left3_algn and not mapped_right5_algn:
        right_overlap = True
    elif not mapped_left3_algn and mapped_right5_algn:
        right_overlap = False
    elif mapped_left3_algn and not mapped_right5_algn:
        right_overlap = False
    else:
        right_overlap = True
        right_overlap &= left3_algn["chrom"] == right5_algn["chrom"]
        right_overlap &= left3_algn["strand"] != right5_algn["strand"]

    same_pair = True
    same_pair &= abs(left5_algn["pos3"] - right3_algn["pos5"]) <= dedup_max_mismatch
    same_pair &= abs(left3_algn["pos5"] - right5_algn["pos3"]) <= dedup_max_mismatch

    if left_overlap & right_overlap & same_pair:
        return 1
    else:
        return 0


def format_pair(
    hic_algn1,
    hic_algn2,
    pair_index,
    report_position="outer",
    report_orientation="pair",
    algn1_pos5=None,
    algn1_pos3=None,
    algn2_pos5=None,
    algn2_pos3=None,
):
    """
    Return a triplet: pair of formatted alignments and pair_index in a walk

    :param hic_algn1: Left alignment forming a pair
    :param hic_algn2: Right alignment forming a pair
    :param algns1: All left read alignments for formal reporting
    :param algns2: All right read alignments for formal reporting
    :param pair_index: Index of the pair
    :param algn1_pos5: Replace reported 5'-position of the alignment 1 with this value
    :param algn1_pos3: Replace reported 3'-position of the alignment 1 with this value
    :param algn2_pos5: Replace reported 5'-position of the alignment 2 with this value
    :param algn2_pos3: Replace reported 3'-position of the alignment 2 with this value

    """
    # Make sure the original data is not modified:
    hic_algn1, hic_algn2 = dict(hic_algn1), dict(hic_algn2)

    # Adjust the 5' and 3'-ends:
    hic_algn1["pos5"] = algn1_pos5 if not algn1_pos5 is None else hic_algn1["pos5"]
    hic_algn1["pos3"] = algn1_pos3 if not algn1_pos3 is None else hic_algn1["pos3"]
    hic_algn2["pos5"] = algn2_pos5 if not algn2_pos5 is None else hic_algn2["pos5"]
    hic_algn2["pos3"] = algn2_pos3 if not algn2_pos3 is None else hic_algn2["pos3"]

    hic_algn1["type"] = (
        "N"
        if not hic_algn1["is_mapped"]
        else "M"
        if not hic_algn1["is_unique"]
        else "U"
    )

    hic_algn2["type"] = (
        "N"
        if not hic_algn2["is_mapped"]
        else "M"
        if not hic_algn2["is_unique"]
        else "U"
    )

    # Change orientation and positioning of pair for reporting:
    # AVAILABLE_REPORT_POSITION    = ["outer", "pair", "read", "walk"]
    # AVAILABLE_REPORT_ORIENTATION = ["pair", "pair", "read", "walk"]
    pair_type = pair_index[1]

    if report_orientation == "read":
        pass
    elif report_orientation == "walk":
        if pair_type == "R2":
            hic_algn1 = flip_orientation(hic_algn1)
            hic_algn2 = flip_orientation(hic_algn2)
        elif pair_type == "R1-2":
            hic_algn2 = flip_orientation(hic_algn2)
    elif report_orientation == "pair":
        if pair_type == "R1" or pair_type == "R1&R2":
            hic_algn2 = flip_orientation(hic_algn2)
        elif pair_type == "R2":
            hic_algn1 = flip_orientation(hic_algn1)
    elif report_orientation == "junction":
        if pair_type == "R1" or pair_type == "R1&R2":
            hic_algn1 = flip_orientation(hic_algn1)
        elif pair_type == "R2":
            hic_algn2 = flip_orientation(hic_algn2)
        else:
            hic_algn1 = flip_orientation(hic_algn1)
            hic_algn2 = flip_orientation(hic_algn2)

    if report_position == "read":
        pass
    elif report_position == "walk":
        if pair_type == "R2":
            hic_algn1 = flip_position(hic_algn1)
            hic_algn2 = flip_position(hic_algn2)
        elif pair_type == "R1-2":
            hic_algn2 = flip_position(hic_algn2)
    elif report_position == "outer":
        if pair_type == "R1" or pair_type == "R1&R2":
            hic_algn2 = flip_position(hic_algn2)
        elif pair_type == "R2":
            hic_algn1 = flip_position(hic_algn1)
    elif report_position == "junction":
        if pair_type == "R1" or pair_type == "R1&R2":
            hic_algn1 = flip_position(hic_algn1)
        elif pair_type == "R2":
            hic_algn2 = flip_position(hic_algn2)
        else:
            hic_algn1 = flip_position(hic_algn1)
            hic_algn2 = flip_position(hic_algn2)

    return [hic_algn1, hic_algn2, pair_index]


def check_pair_order(algn1, algn2, chrom_enum):
    """
    Check if a pair of alignments has the upper-triangular order or
    has to be flipped.
    """

    # First, the pair is flipped according to the type of mapping on its sides.
    # Later, we will check it is mapped on both sides and, if so, flip the sides
    # according to these coordinates.
    has_correct_order = (algn1["is_mapped"], algn1["is_unique"]) <= (
        algn2["is_mapped"],
        algn2["is_unique"],
    )

    # If a pair has coordinates on both sides, it must be flipped according to
    # its genomic coordinates.
    if (algn1["chrom"] != pairsam_format.UNMAPPED_CHROM) and (
        algn2["chrom"] != pairsam_format.UNMAPPED_CHROM
    ):
        has_correct_order = (chrom_enum[algn1["chrom"]], algn1["pos"]) <= (
            chrom_enum[algn2["chrom"]],
            algn2["pos"],
        )

    return has_correct_order


######################
### Output utilities:
######################


def write_all_algnments(readID, all_algns1, all_algns2, out_file):
    """
    Debug utility that outputs all alignments in .bam file before parsing walks/pairs
    """
    for side_idx, all_algns in enumerate((all_algns1, all_algns2)):
        out_file.write(readID)
        out_file.write("\t")
        out_file.write(str(side_idx + 1))
        out_file.write("\t")
        for algn in sorted(all_algns, key=lambda x: x["dist_to_5"]):
            out_file.write(algn["chrom"])
            out_file.write("\t")
            out_file.write(str(algn["pos5"]))
            out_file.write("\t")
            out_file.write(algn["strand"])
            out_file.write("\t")
            out_file.write(str(algn["mapq"]))
            out_file.write("\t")
            out_file.write(str(algn["cigar"]))
            out_file.write("\t")
            out_file.write(str(algn["dist_to_5"]))
            out_file.write("\t")
            out_file.write(str(algn["dist_to_5"] + algn["algn_read_span"]))
            out_file.write("\t")
            out_file.write(str(algn["matched_bp"]))
            out_file.write("\t")

        out_file.write("\n")


def write_pairsam(
    algn1,
    algn2,
    readID,
    pair_index,
    sams1,
    sams2,
    out_file,
    drop_readid,
    drop_seq,
    drop_sam,
    add_pair_index,
    add_columns,
):
    """
    Write output pairsam.
    Note: SAM is already tab-separated and
    any printable character between ! and ~ may appear in the PHRED field!
    (http://www.ascii-code.com/)
    Thus, use the vertical tab character to separate fields!
    """
    cols = [
        "." if drop_readid else readID,
        algn1["chrom"],
        str(algn1["pos"]),
        algn2["chrom"],
        str(algn2["pos"]),
        algn1["strand"],
        algn2["strand"],
        algn1["type"] + algn2["type"],
    ]

    if not drop_sam:
        for sams in [sams1, sams2]:
            if drop_seq:
                for sam in sams:
                    sam.query_qualities = ""
                    sam.query_sequence = ""
            cols.append(
                pairsam_format.INTER_SAM_SEP.join(
                    [
                        sam.to_string().replace(
                            "\t", pairsam_format.SAM_SEP
                        )  # String representation of pysam alignment
                        + pairsam_format.SAM_SEP
                        + "Yt:Z:"
                        + algn1["type"]
                        + algn2["type"]
                        for sam in sams
                    ]
                )
            )

    if add_pair_index:
        cols.append(str(pair_index[0]))
        cols.append(pair_index[1])

    for col in add_columns:
        # use get b/c empty alignments would not have sam tags (NM, AS, etc)
        cols.append(str(algn1.get(col, "")))
        cols.append(str(algn2.get(col, "")))

    out_file.write(pairsam_format.PAIRSAM_SEP.join(cols) + "\n")
