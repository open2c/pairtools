"""
Set of functions used for pairsam parse, migrated from pairtools/pairtools_parse.py

Parse operates with several basic data types:

I. pysam-based:
    1. **sam entry** is a continuous aligned fragment of the read mapped to certain location in the genome.
        Because we read sam entries from .sam/.bam files automatically with modified pysam,
        each sam entry is in fact special AlignedSegmentPairtoolized Cython object
        that has alignment attributes and can be easily accessed from Python.
    
        Sam entries are gathered into reads by `push_pysam` function.
    
    2. **read** is a collection of sam entries corresponding to a single Hi-C molecule.
        It is represented by three variables:
        readID, sams1 and sams2, which keep left and right sam entries, correspondingly.
        Read is populated from the stream of sam entries on a fly, the process happenning
        in `streaming_classify` function.

II. python-based data types are parsed from pysam-based ones:

    1. **alignment** is a continuous aligned fragment represented as dictionary with relevant fields,
        such as "chrom", "pos5", "pos3", "strand", "type", etc.
        
        `empty_alignment` creates empty alignment,
        `parse_pysam_entry` create new alignmetns from pysam entries,
        `mask_alignment` clears some fields of the alignment to match the default "unmapped" state.
        
        `flip_alignment`, `flip_orientation` and `flip_ends` are useful functions that help to orient alignments.

    2. **pair** of two alignments is represented by three variables:
        algn1 (left alignment), algn2 (right alignment) and junction_index.
        Pairs are obtained by `parse_read` or `parse2_read`.
        Additionally, these functions also output all alignments for each side.

"""

from . import _pairsam_format

def streaming_classify(
    instream,
    outstream,
    chromosomes,
    out_alignments_stream,
    out_stat,
    **kwargs
):
    """
    Parse input sam file into individual reads, pairs, walks,
    then write to the outstream(s).

    Additional kwargs:
        min_mapq,
        drop_readid,
        drop_seq,
        drop_sam,
        add_junction_index,
        add_columns,
        report_alignment_end,
        max_inter_align_gap
    parse:
        max_molecule_size,
        walks_policy
    parse2:
        max_fragment_size,
        single_end,
        report_position,
        report_orientation

    """

    parse2 = kwargs.get("parse2", False)

    ### Store output parameters in a usable form:
    chrom_enum = dict(
        zip(
            [_pairsam_format.UNMAPPED_CHROM] + list(chromosomes),
            range(len(chromosomes) + 1),
        )
    )
    add_columns = kwargs.get("add_columns", [])
    sam_tags = [col for col in add_columns if len(col) == 2 and col.isupper()]
    store_seq = "seq" in add_columns

    ### Compile readID transformation:
    readID_transform = kwargs.get("readid_transform", None)
    if readID_transform is not None:
        readID_transform = compile(readID_transform, "<string>", "eval")

    ### Prepare for iterative parsing of the input stream
    # Each read is represented by readID, sams1 (left alignments) and sams2 (right alignments)
    readID = "" # Read id of the current read
    sams1 = []  # Placeholder for the left alignments
    sams2 = []  # Placeholder for the right alignments
    # Each read is comprised of multiple alignments, or sam entries:
    sam_entry = ""  # Placeholder for each aligned segment
    # Keep the id of the previous sam entry to detect when the read is completely populated:
    prev_readID = ""  # Placeholder for the read id

    ### Iterate over input pysam:
    instream = iter(instream)
    while sam_entry is not None:
        sam_entry = next(instream, None)

        readID = sam_entry.query_name if sam_entry else None
        if readID_transform is not None and readID is not None:
            readID = eval(readID_transform)

        # Read is fully populated, then parse and write:
        if not (sam_entry) or ((readID != prev_readID) and prev_readID):

            ### Parse
            if not parse2: # regular parser:
                pairstream = parse_read(
                    sams1,
                    sams2,
                    kwargs["min_mapq"],
                    kwargs["max_molecule_size"],
                    kwargs["max_inter_align_gap"],
                    kwargs["walks_policy"],
                    sam_tags,
                    store_seq
                )
            else:  # parse2 parser:
                pairstream = parse2_read(
                    sams1,
                    sams2,
                    kwargs["min_mapq"],
                    kwargs["max_inter_align_gap"],
                    kwargs["max_fragment_size"],
                    kwargs["single_end"],
                    kwargs["report_position"],
                    kwargs["report_orientation"],
                    sam_tags,
                    store_seq
                )

            ### Write:
            read_has_alignments = False
            for (
                algn1,
                algn2,
                all_algns1,
                all_algns2,
                junction_index,
            ) in pairstream:
                read_has_alignments = True

                if kwargs["report_alignment_end"] == "5":
                    algn1["pos"] = algn1["pos5"]
                    algn2["pos"] = algn2["pos5"]
                else:
                    algn1["pos"] = algn1["pos3"]
                    algn2["pos"] = algn2["pos3"]

                if not kwargs["no_flip"]:
                    flip_pair = not check_pair_order(algn1, algn2, chrom_enum)
                    if flip_pair:
                        algn1, algn2 = algn2, algn1
                        sams1, sams2 = sams2, sams1

                write_pairsam(
                    algn1,
                    algn2,
                    prev_readID,
                    junction_index,
                    sams1,
                    sams2,
                    outstream,
                    kwargs["drop_readid"],
                    kwargs["drop_seq"],
                    kwargs["drop_sam"],
                    kwargs["add_junction_index"],
                    kwargs["add_columns"]
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
                    prev_readID, all_algns1, all_algns2, out_alignments_stream
                )

            # Empty read after writing:
            sams1.clear()
            sams2.clear()

        if sam_entry is not None:
            push_pysam(sam_entry, sams1, sams2)
            prev_readID = readID


####################
### Pysam utilities:
####################

def push_pysam(sam_entry, sams1, sams2):
    """Parse pysam AlignedSegment (sam) into pairtools sams entry"""
    flag = sam_entry.flag
    if (flag & 0x40) != 0:
        sams1.append(sam_entry)  # Forward read, or first read in a pair
    else:
        sams2.append(sam_entry)  # Reverse read, or mate pair
    return


############################
### Alignment utilities: ###
############################

def empty_alignment():
    return {
        "chrom": _pairsam_format.UNMAPPED_CHROM,
        "pos5": _pairsam_format.UNMAPPED_POS,
        "pos3": _pairsam_format.UNMAPPED_POS,
        "pos": _pairsam_format.UNMAPPED_POS,
        "strand": _pairsam_format.UNMAPPED_STRAND,
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
    }

def parse_pysam_entry(
    sam, min_mapq, sam_tags=None, store_seq=False, report_3_alignment_end=False
):
    """Parse alignments from pysam AlignedSegment entry
    :param sam: input pysam AlignedSegment entry
    :param min_mapq: minimal MAPQ to consider as a proper alignment
    :param sam_tags: list of sam tags to store
    :param store_seq: if True, the sequence will be parsed and stored in the output
    :param report_3_alignment_end: if True, 3'-end of alignment will be reported as position (will be deprecated)
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
                # print(cigar['algn_ref_span'])
                # Note that pysam output is zero-based, thus add +1:
                pos5 = sam.reference_start + 1
                pos3 = sam.reference_start + cigar["algn_ref_span"]
            else:
                pos5 = sam.reference_start + cigar["algn_ref_span"]
                # Note that pysam output is zero-based, thus add +1:
                pos3 = sam.reference_start + 1
            # print(pos5, pos3)

        else:
            chrom = _pairsam_format.UNMAPPED_CHROM
            strand = _pairsam_format.UNMAPPED_STRAND
            pos5 = _pairsam_format.UNMAPPED_POS
            pos3 = _pairsam_format.UNMAPPED_POS
    else:
        chrom = _pairsam_format.UNMAPPED_CHROM
        strand = _pairsam_format.UNMAPPED_STRAND
        pos5 = _pairsam_format.UNMAPPED_POS
        pos3 = _pairsam_format.UNMAPPED_POS

        dist_to_5 = 0
        dist_to_3 = 0

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
    algn["chrom"] = _pairsam_format.UNMAPPED_CHROM
    algn["pos5"] = _pairsam_format.UNMAPPED_POS
    algn["pos3"] = _pairsam_format.UNMAPPED_POS
    algn["pos"] = _pairsam_format.UNMAPPED_POS
    algn["strand"] = _pairsam_format.UNMAPPED_STRAND

    return algn

def flip_alignment(hic_algn):
    """
    Flip a single alignment as if it was sequenced from the opposite end
    :param hic_algn: Alignment to be modified
    :return:
    """
    hic_algn = dict(hic_algn)  # overwrite the variable with the copy of dictionary
    hic_algn["pos5"], hic_algn["pos3"] = hic_algn["pos3"], hic_algn["pos5"]
    hic_algn["strand"] = "+" if hic_algn["strand"] == "-" else "-"
    return hic_algn

def flip_orientation(hic_algn):
    """
    Flip orientation of a single alignment
    :param hic_algn: Alignment to be modified
    :return:
    """
    hic_algn = dict(hic_algn)  # overwrite the variable with the copy of dictionary
    hic_algn["strand"] = "+" if hic_algn["strand"] == "-" else "-"
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
    store_seq
):
    """
    Parse sam entries corresponding to a single read (or Hi-C molecule)
    into pairs of alignments.

    Returns
    -------
    algn1, algn2: dict
        Two alignments selected for reporting as a Hi-C pair.
    algns1, algns2
        All alignments, sorted according to their order in on a read.
    junction_index
        Junction index of a pair in the molecule.
    """

    # Check if there is at least one sam entry per side:
    if (len(sams1) == 0) or (len(sams2) == 0):
        algns1 = [empty_alignment()]
        algns2 = [empty_alignment()]
        algns1[0]["type"] = "X"
        algns2[0]["type"] = "X"
        junction_index = "1u"
        return [[algns1[0], algns2[0], algns1, algns2, junction_index]]

    # Generate a sorted, gap-filled list of all alignments
    algns1 = [ parse_pysam_entry(sam, min_mapq, sam_tags, store_seq) for sam in sams1 ]
    algns2 = [ parse_pysam_entry(sam, min_mapq, sam_tags, store_seq) for sam in sams2 ]

    algns1 = sorted(algns1, key=lambda algn: algn["dist_to_5"])
    algns2 = sorted(algns2, key=lambda algn: algn["dist_to_5"])

    if max_inter_align_gap is not None:
        _convert_gaps_into_alignments(algns1, max_inter_align_gap)
        _convert_gaps_into_alignments(algns2, max_inter_align_gap)

    # By default, assume each molecule is a single pair with single unconfirmed junction:
    hic_algn1 = algns1[0]
    hic_algn2 = algns2[0]
    junction_index = "1u"

    # Define the type of alignment on each side:
    is_chimeric_1 = len(algns1) > 1
    is_chimeric_2 = len(algns2) > 1

    # Parse chimeras
    if is_chimeric_1 or is_chimeric_2:

        # Report all the linear alignments in a read pair
        if walks_policy == "all":
            # Report linear alignments after deduplication of complex walks with default settings:
            return parse_complex_walk(algns1, algns2, max_molecule_size,
                                      report_position="outer",
                                      report_orientation="pair")

        elif walks_policy in ['mask', '5any', '5unique', '3any', '3unique']:
            # Report only two alignments for a read pair
            rescued_linear_side = rescue_walk(algns1, algns2, max_molecule_size)

            # Walk was rescued as a simple walk:
            if rescued_linear_side is not None:
                junction_index = f'1{"f" if rescued_linear_side==1 else "r"}'
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

                # lower-case reported walks on the chimeric side
                if walks_policy != "mask":
                    if is_chimeric_1:
                        hic_algn1 = dict(hic_algn1)
                        hic_algn1["type"] = hic_algn1["type"].lower()
                    if is_chimeric_2:
                        hic_algn2 = dict(hic_algn2)
                        hic_algn2["type"] = hic_algn2["type"].lower()

        else:
            raise ValueError(f"Walks policy {walks_policy} is not supported.")

    return [[hic_algn1, hic_algn2, algns1, algns2, junction_index]]


#### parse2 parser:
def parse2_read(
    sams1,
    sams2,
    min_mapq,
    max_inter_align_gap,
    max_fragment_size,
    single_end,
    report_position="outer",
    report_orientation="pair",
    sam_tags=[],
    store_seq=False
):
    """
    Parse sam entries corresponding to a Hi-C molecule into alignments in parse2 mode
    for a Hi-C pair.
    Returns
    -------
    algn1, algn2: dict
        Two alignments selected for reporting as a Hi-C pair.
    algns1, algns2
        All alignments, sorted according to their order in on a read.
    junction_index
        Junction index of a pair in the molecule.
    """

    # Single-end mode:
    if single_end:
        # Generate a sorted, gap-filled list of all alignments
        algns1 = [parse_pysam_entry(sam, min_mapq, sam_tags, store_seq) for sam in sams1]
        algns1 = sorted(algns1, key=lambda algn: algn["dist_to_5"])
        if max_inter_align_gap is not None:
            _convert_gaps_into_alignments(algns1, max_inter_align_gap)

        algns2 = [empty_alignment()]  # Empty alignment dummy

        if len(algns1) > 1:
            # Look for ligation junction, and report linear alignments after deduplication of complex walks:
            # (Note that coordinate system for single-end reads does not change the behavior)
            return parse_complex_walk(
                algns1, algns2, max_fragment_size, report_position, report_orientation
            )
        else:
            # If no additional information, we assume each molecule is a single ligation with single unconfirmed junction:
            algn2 = algns2[0]
            if report_orientation == "walk":
               algn2 = flip_orientation(algn2)
            if report_position == "walk":
                algn2 = flip_position(algn2)
            return [[algns1[0], algn2, algns1, algns2, "1u"]]

    # Paired-end mode:
    else:
        # Check if there is at least one SAM entry per side:
        if (len(sams1) == 0) or (len(sams2) == 0):
            algns1 = [empty_alignment()]
            algns2 = [empty_alignment()]
            algns1[0]["type"] = "X"
            algns2[0]["type"] = "X"
            return [[algns1[0], algns2[0], algns1, algns2, "1u"]]

        # Generate a sorted, gap-filled list of all alignments
        algns1 = [parse_pysam_entry(sam, min_mapq, sam_tags, store_seq) for sam in sams1]
        algns2 = [parse_pysam_entry(sam, min_mapq, sam_tags, store_seq) for sam in sams2]

        algns1 = sorted(algns1, key=lambda algn: algn["dist_to_5"])
        algns2 = sorted(algns2, key=lambda algn: algn["dist_to_5"])

        if max_inter_align_gap is not None:
            _convert_gaps_into_alignments(algns1, max_inter_align_gap)
            _convert_gaps_into_alignments(algns2, max_inter_align_gap)

        is_chimeric_1 = len(algns1) > 1
        is_chimeric_2 = len(algns2) > 1

        if is_chimeric_1 or is_chimeric_2:
            # If at least one side is chimera, we must look for ligation junction, and
            # report linear alignments after deduplication of complex walks:
            return parse_complex_walk(
                algns1, algns2, max_fragment_size, report_position, report_orientation
            )
        else:
            # If no additional information, we assume each molecule is a single ligation with single unconfirmed junction:
            algn2 = algns2[0]
            if report_orientation == "walk":
               algn2 = flip_orientation(algn2)
            if report_position == "walk":
                algn2 = flip_position(algn2)
            return [[algns1[0], algn2, algns1, algns2, "1u"]]


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
    a) the 3'-end alignment on one side maps to the same chromosome as the
    alignment fully covering the other side (i.e. the linear alignment)
    b) the two alignments point towards each other on the chromosome
    c) the distance between the outer ends of the two alignments is below
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
        if first_read_is_chimeric:
            # changing the type of the 3' alignment on side 1, does not show up
            # in the output
            algns1[1]["type"] = "X"
            algns2[0]["type"] = "R"
            return 1
        else:
            algns1[0]["type"] = "R"
            # changing the type of the 3' alignment on side 2, does not show up
            # in the output
            algns2[1]["type"] = "X"
            return 2
    else:
        return None

def _convert_gaps_into_alignments(sorted_algns, max_inter_align_gap):
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


#### Complex walks parser:
def parse_complex_walk(
    algns1,
    algns2,
    max_fragment_size,
    report_position,
    report_orientation,
    allowed_offset=3
):
    """
    Parse a set of ligations that appear as a complex walk.

    If the reads are long enough, the reverse read might read through the forward read's meaningful part.
    And if one of the reads contains ligation junction, this might lead to reporting a fake contact!
    Thus, the pairs of contacts that overlap between forward and reverse reads are paired-end duplicates.
    This complex walk parser treats these cases and reports only unique pairs of alignments as contacts.

    :param algns1: List of sequential forwards alignments
    :param algns2: List of sequential reverse alignments
    :param max_fragment_size:
    :param report_position:
    :param report_orientation:
    :param allowed_offset: the number of basepairs that are allowed at the ends of alignments to detect overlaps

    :return: list of all the pairs after paired-end deduplication.

    Illustration of the algorithm inner working.

     Forward read:                            Reverse read:
    ---------------------->       <-----------------------
             algns1                        algns2
    5---3_5---3_5---3_5---3        3---5_3---5_3---5_3---5
        fIII  fII   fI                 rI    rII  rIII
          junctions                       junctions

    Alignment is a bwa mem reported hit. After parsing of bam file, all the alignments are reported in
    sequential order as algns1 for forward and algns2 for reverse reads.
    Junction is a sequential pair of linear alignments reported as chimera at forward or reverse read.

    Let's consider the case if n_algns1 >= 2 on forward read and n_algns2 >= 2 on reverse read.
    We start looking for overlapping pairs of linear alignments from the ends of reads.

    The procedure of iterative search of overlap:
      1. Take the last 3' junction on the forward read (fI, or current_forward_junction)
          and the last 3' junction on reverse read (rI, or current_reverse_junction).
      2. Compare fI and rI (pairs_do_overlap).
          If successful, we found the overlap, add it to the output list.
          If not successful, go to p.3.
      3. Take the next pair of linear alignments of reverse read (rII), i.e. shift current_reverse_junction by one.
      4. Check that this pair can form a potential overlap with fI:
            the number of junctions downstream from fI on forward read should not be less than
            the number of junctions upstream from rII on reverse read.
         If the potential overlap can be formed, go to p. 5.
         If it cannot be formed, no other overlap in this complex walk is possible. Exit.
      5. Compare the current pair of junctions on forward and reverse reads.
         If comparison fails, go to p. 3, i.e. take the next pair of linear alignments of reverse read (rIII).
         If comparison is successful, check that junctions downstream from fI overlap with the junctions upstream from rII.
             If yes, add them all to the output list.
             If not, we do not have an overlap, repeat p. 3.

    Note that we do not need to shift forward read, because biologically overlap can only happen
    when both ends of forward and reverse read are involved, and shifting one of them is enough.
    """

    AVAILABLE_REPORT_POSITION    = ["outer", "junction", "read", "walk"]
    assert report_position in AVAILABLE_REPORT_POSITION, (
        f"Cannot report position {report_position}, as it is not implemented"
        f'Available choices are: {", ".join(AVAILABLE_REPORT_POSITION)}'
    )

    AVAILABLE_REPORT_ORIENTATION = ["pair", "junction", "read", "walk"]
    assert report_orientation in AVAILABLE_REPORT_ORIENTATION, (
        f"Cannot report orientation {report_orientation}, as it is not implemented"
        f'Available choices are: {", ".join(AVAILABLE_REPORT_ORIENTATION)}'
    )

    n_algns1 = len(algns1)
    n_algns2 = len(algns2)

    ### Complex walk parser algorithm ###

    # Storage for the final contacts:
    final_contacts = []

    # Initialize some useful variables:
    current_forward_junction = current_reverse_junction = 1  # p. 1, initialization
    remaining_forward_junctions = n_algns1 - 1 # Number of possible junctions remaining on forward read
    remaining_reverse_junctions = n_algns2 - 1 # Number of possible junctions remaining on reverse read
    checked_reverse_junctions = 0  # Number of checked junctions on reverse read (from the end of read)
    is_overlap = False

    # Iterative search of overlaps between forward and reverse alignments.
    # If both sides have more than 2 alignments, then check if there are overlapping forward and reverse alignments pairs:
    if (n_algns1 >= 2) and (n_algns2 >= 2):

        # Loop through all alignment pairs and check for overlaps:
        while (remaining_forward_junctions > checked_reverse_junctions) and (remaining_reverse_junctions > 0):

            # Check if current pairs of junctions overlap:
            is_overlap = pairs_do_overlap(
                (
                    algns1[-current_forward_junction - 1],
                    algns1[-current_forward_junction],
                ),
                (
                    algns2[-current_reverse_junction - 1],
                    algns2[-current_reverse_junction],
                ),
                allowed_offset,
            )

            # There is a potential overlap, we need to check whether it's consistent,
            # i.e. that the remaining pairs of forward downstream and reverse upstream junctions overlap as well:
            if is_overlap:
                last_idx_forward_temp = current_forward_junction
                last_idx_reverse_temp = current_reverse_junction
                checked_reverse_temp = checked_reverse_junctions
                # loop over all forward downstream and reverse upstream junctions:
                while is_overlap and (checked_reverse_temp > 0):
                    last_idx_forward_temp += 1
                    last_idx_reverse_temp -= 1
                    is_overlap &= pairs_do_overlap(
                        (
                            algns1[-last_idx_forward_temp - 1],
                            algns1[-last_idx_forward_temp],
                        ),
                        (
                            algns2[-last_idx_reverse_temp - 1],
                            algns2[-last_idx_reverse_temp],
                        ),
                        allowed_offset,
                    )
                    checked_reverse_temp -= 1
                # all the checks have passed, no need to check for another hit:
                if is_overlap:
                    current_reverse_junction += 1
                    break

            # p. 3: shift the reverse junction pointer by one
            current_reverse_junction += 1
            checked_reverse_junctions += 1
            remaining_reverse_junctions -= 1

        # No overlap found, roll the current_idx_reverse back to the initial value:
        if not is_overlap:
            current_reverse_junction = 1

    # If there are less than 2 chimeras in either forward or reverse read, or no overlapping junctions found,
    # then current_reverse_junction is 1, and we check whether the last alignments of forward and reverse reads overlap.
    if current_reverse_junction == 1:
        last_reported_alignment_forward = last_reported_alignment_reverse = 1
        # If the last alignments on forward and reverse overlap, then report the last pairs of junctions on each side:
        if ends_do_overlap(algns1[-1], algns2[-1], max_fragment_size, allowed_offset):

            # Report the last of multiple alignments on forward read and single alignment on reverse:
            if (n_algns1 >= 2):
                push_pair(
                    final_contacts,
                    algns1[-2],
                    algns1[-1],
                    algns1,
                    algns2,
                    junction_index=f"{len(algns1)-1}f",
                    algn2_pos3=algns2[-1]["pos5"],
                    report_position=report_position,
                    report_orientation=report_orientation
                )
                last_reported_alignment_forward = 2

            # Single alignment on forward read and multiple alignments on reverse:
            if (n_algns2 >= 2):
                push_pair(
                    final_contacts,
                    algns2[-1],
                    algns2[-2],
                    algns1,
                    algns2,
                    junction_index=f"{len(algns1)}r",
                    algn1_pos3=algns1[-1]["pos5"],
                    report_position=report_position,
                    report_orientation=report_orientation
                )
                last_reported_alignment_reverse = 2

            # Note that if n_algns1==n_algns2==1 and alignments overlap, then we don't need to check,
            # it's a non-ligated DNA fragment that we don't report.

        # If end alignments do not overlap, then there is no evidence of ligation junction for the pair,
        # report a regular pair:
        else:
            push_pair(
                final_contacts,
                algns1[-1],
                algns2[-1],
                algns1,
                algns2,
                junction_index=f"{len(algns1)}u",
                report_position=report_position,
                report_orientation=report_orientation
            )

    # If we have an overlap of junctions:
    else:
        last_reported_alignment_forward = last_reported_alignment_reverse = current_reverse_junction

    # Report all unique alignments on forward read (sequential):
    for i in range(0, n_algns1 - last_reported_alignment_forward):
        push_pair(
            final_contacts,
            algns1[i],
            algns1[i + 1],
            algns1,
            algns2,
            junction_index=f"{i + 1}f",
            report_position=report_position,
            report_orientation=report_orientation
        )

    # Report the pairs where both forward alignments overlap reverse:
    for i_overlapping in range(current_reverse_junction - 1):
        idx_forward = n_algns1 - current_reverse_junction + i_overlapping
        idx_reverse = n_algns2 - 1 - i_overlapping
        push_pair(
            final_contacts,
            algns1[idx_forward],
            algns1[idx_forward + 1],
            algns1,
            algns2,
            junction_index=f"{idx_forward + 1}b",
            algn2_pos3=algns2[idx_reverse - 1]["pos5"],
            report_position=report_position,
            report_orientation=report_orientation
        )

    # Report all the sequential chimeric pairs in the reverse read, but not the overlap:
    if report_position in ['walk']: # Report from 3'-end to 5'-end of reverse read to keep the walk order
        for i in range(0, min(current_reverse_junction, n_algns2 - last_reported_alignment_reverse))[::-1]:
            # Determine the junction index depending on what is the overlap:
            if current_reverse_junction > 1:
                junction_index = n_algns1 + min(current_reverse_junction,
                                                n_algns2 - last_reported_alignment_reverse) - i - 1
            else:
                junction_index = n_algns1 + min(current_reverse_junction,
                                                n_algns2 - last_reported_alignment_reverse) - i

            push_pair(
                final_contacts,
                algns2[i+1],
                algns2[i],
                algns1,
                algns2,
                junction_index=f"{junction_index}r",
                report_position=report_position,
                report_orientation=report_orientation
            )

    else: # Report from 5'-end to 3'-end of reverse read to keep the read order
        for i in range(0, min(current_reverse_junction, n_algns2 - last_reported_alignment_reverse)):
            # Determine the junction index depending on what is the overlap:
            if current_reverse_junction > 1:
                junction_index = n_algns1 + min(current_reverse_junction,
                                                n_algns2 - last_reported_alignment_reverse) - i - 1
            else:
                junction_index = n_algns1 + min(current_reverse_junction,
                                                n_algns2 - last_reported_alignment_reverse) - i

            push_pair(
                final_contacts,
                algns2[i+1],
                algns2[i],
                algns1,
                algns2,
                junction_index=f"{junction_index}r",
                report_position=report_position,
                report_orientation=report_orientation
            )

    # Sort the pairs according to the order of appearance in the reads.
    # Take the junction index (last element in each entry from its end),
    # and put forward reads first, then the reverse reads:
    final_contacts.sort(key=lambda x: int(x[-1][:-1]))
    return final_contacts


### Additional functions for complex walks rescue ###
def ends_do_overlap(algn1, algn2, max_fragment_size=500, allowed_offset=5):
    """
    Two ends of alignments overlap if:
     1) they are from the same chromosome,
     2) map in the opposite directions,
     3) the distance between the outer ends of the two alignments is below the specified max_fragment_size,
     4) the distance between the outer ends of the two alignments is above the maximum alignment size.
    (4) guarantees that the alignments point towards each other on the chromosomes.

    Allowed offset is for the cases when few nucleotides are mismapped by bwa at the ends of chimeric parts.

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

    do_overlap &= distance_outer_ends <= max_fragment_size + allowed_offset
    do_overlap &= distance_outer_ends >= min_algn_size - allowed_offset

    if do_overlap:
        return 1
    return 0


def pairs_do_overlap(algns1, algns2, allowed_offset=5):
    """
    Forward read:                             Reverse read:
    ----------------------->      <------------------------
             algns1                        algns2
    5----------3_5----------3     3----------5_3----------5
    algn1_chim5   algn1_chim3     algn2_chim3   algn2_chim5
    chim_left     chim_right      chim_left     chim_right

    Two pairs of alignments overlap if:
    1) algn1_chim5 and algn2_chim3 originate from the same region (chim_left),
    2) algn1_chim3 and algn2_chim5 originate from the same region (chim_right).
    or:
    3) pos3 of algn1_chim5 is close to pos3 of algn2_chim3,
    4) pos5 of algn1_chim3 is close to pos5 of algn2_chim5.

    Return: 1 of the pairs of alignments are overlaps,
            0 if they are not.
    """

    # Some assignments to simplify the code
    algn1_chim5 = algns1[0]
    algn1_chim3 = algns1[1]
    algn2_chim5 = algns2[0]
    algn2_chim3 = algns2[1]

    # We assume that successful alignment cannot be an overlap with unmapped or multi-mapped region
    mapped_algn1_chim5 = algn1_chim5["is_mapped"] and algn1_chim5["is_unique"]
    mapped_algn1_chim3 = algn1_chim3["is_mapped"] and algn1_chim3["is_unique"]
    mapped_algn2_chim5 = algn2_chim5["is_mapped"] and algn2_chim5["is_unique"]
    mapped_algn2_chim3 = algn2_chim3["is_mapped"] and algn2_chim3["is_unique"]

    if not mapped_algn1_chim5 and not mapped_algn2_chim3:
        chim_left_overlap = True
    elif not mapped_algn1_chim5 and mapped_algn2_chim3:
        chim_left_overlap = False
    elif mapped_algn1_chim5 and not mapped_algn2_chim3:
        chim_left_overlap = False
    else:
        chim_left_overlap = True
        chim_left_overlap &= algn1_chim5["chrom"] == algn2_chim3["chrom"]
        chim_left_overlap &= algn1_chim5["strand"] != algn2_chim3["strand"]

    if not mapped_algn1_chim3 and not mapped_algn2_chim5:
        chim_right_overlap = True
    elif not mapped_algn1_chim3 and mapped_algn2_chim5:
        chim_right_overlap = False
    elif mapped_algn1_chim3 and not mapped_algn2_chim5:
        chim_right_overlap = False
    else:
        chim_right_overlap = True
        chim_right_overlap &= algn1_chim3["chrom"] == algn2_chim5["chrom"]
        chim_right_overlap &= algn1_chim3["strand"] != algn2_chim5["strand"]

    same_junction = True
    same_junction &= abs(algn1_chim5["pos3"] - algn2_chim3["pos5"]) <= allowed_offset
    same_junction &= abs(algn1_chim3["pos5"] - algn2_chim5["pos3"]) <= allowed_offset

    if chim_left_overlap & chim_right_overlap & same_junction:
        return 1
    else:
        return 0


def push_pair(
    final_contacts,
    hic_algn1,
    hic_algn2,
    algns1,
    algns2,
    junction_index,
    report_position="outer",
    report_orientation="pair",
    algn1_pos5=None,
    algn1_pos3=None,
    algn2_pos5=None,
    algn2_pos3=None,
):
    """
    Push a pair of alignments into final list of contacts.

    :param final_contacts: List that will be updated

    :param hic_algn1: Left alignment forming a pair
    :param hic_algn2: Right alignment forming a pair

    :param algns1: All forward read alignments for formal reporting
    :param algns2: All reverse read alignments for formal reporting

    :param junction_index: Index of the junction

    :param algn1_pos5: Replace reported 5'-position of the alignment 1 with this value
    :param algn1_pos3: Replace reported 3'-position of the alignment 1 with this value
    :param algn2_pos5: Replace reported 5'-position of the alignment 2 with this value
    :param algn2_pos3: Replace reported 3'-position of the alignment 2 with this value

    :return: 0 if successful
    """

    # Overwrite the variables with copies of dictionaries
    # to make sure the original data is not modified:
    hic_algn1, hic_algn2 = dict(hic_algn1), dict(hic_algn2)

    # Adjust the 5' and 3'-ends:
    hic_algn1["pos5"] = algn1_pos5 if not algn1_pos5 is None else hic_algn1["pos5"]
    hic_algn1["pos3"] = algn1_pos3 if not algn1_pos3 is None else hic_algn1["pos3"]
    hic_algn2["pos5"] = algn2_pos5 if not algn2_pos5 is None else hic_algn2["pos5"]
    hic_algn2["pos3"] = algn2_pos3 if not algn2_pos3 is None else hic_algn2["pos3"]

    hic_algn1["type"] = "N" if not hic_algn1["is_mapped"] else \
                        "M" if not hic_algn1["is_unique"] else \
                        "U"

    hic_algn2["type"] = "N" if not hic_algn2["is_mapped"] else \
                        "M" if not hic_algn2["is_unique"] else \
                        "U"

    # Determine orientation and ositioning of the pair:
    # AVAILABLE_REPORT_POSITION    = ["outer", "junction", "read", "walk"]
    # AVAILABLE_REPORT_ORIENTATION = ["pair", "junction", "read", "walk"]
    pair_type = junction_index[-1]

    if report_orientation=="read":
        pass
    elif report_orientation=="walk":
        if pair_type=="r":
            hic_algn1 = flip_orientation(hic_algn1)
            hic_algn2 = flip_orientation(hic_algn2)
        elif pair_type=="u":
            hic_algn2 = flip_orientation(hic_algn2)
    elif report_orientation=="pair":
        if pair_type in ["f", "r"]:
            hic_algn2 = flip_orientation(hic_algn2)
    elif report_orientation=="junction":
        if pair_type in ["f", "r"]:
            hic_algn1 = flip_orientation(hic_algn1)
        else:
            hic_algn1 = flip_orientation(hic_algn1)
            hic_algn2 = flip_orientation(hic_algn2)

    if report_position=="read":
        pass
    elif report_position=="walk":
        if pair_type=="r":
            hic_algn1 = flip_position(hic_algn1)
            hic_algn2 = flip_position(hic_algn2)
        elif pair_type=="u":
            hic_algn2 = flip_position(hic_algn2)
    elif report_position=="outer":
        if pair_type in ["f", "r"]:
            hic_algn2 = flip_position(hic_algn2)
    elif report_position=="junction":
        if pair_type in ["f", "r"]:
            hic_algn1 = flip_position(hic_algn1)
        else:
            hic_algn1 = flip_position(hic_algn1)
            hic_algn2 = flip_position(hic_algn2)

    final_contacts.append([hic_algn1, hic_algn2, algns1, algns2, junction_index])

    return 0


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
    if (algn1["chrom"] != _pairsam_format.UNMAPPED_CHROM) and (
        algn2["chrom"] != _pairsam_format.UNMAPPED_CHROM
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
    junction_index,
    sams1,
    sams2,
    out_file,
    drop_readid,
    drop_seq,
    drop_sam,
    add_junction_index,
    add_columns,
):
    """
    SAM is already tab-separated and
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
                    sam.query_qualities = ''
                    sam.query_sequence = ''
            cols.append(
                _pairsam_format.INTER_SAM_SEP.join(
                    [
                        sam.to_string().replace(
                            "\t", _pairsam_format.SAM_SEP
                        )  # String representation of pysam alignment
                        + _pairsam_format.SAM_SEP
                        + "Yt:Z:"
                        + algn1["type"]
                        + algn2["type"]
                        for sam in sams
                    ]
                )
            )

    if add_junction_index:
        cols.append(junction_index)

    for col in add_columns:
        # use get b/c empty alignments would not have sam tags (NM, AS, etc)
        cols.append(str(algn1.get(col, "")))
        cols.append(str(algn2.get(col, "")))

    out_file.write(_pairsam_format.PAIRSAM_SEP.join(cols) + "\n")
