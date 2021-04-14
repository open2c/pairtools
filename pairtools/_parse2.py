"""
Set of functions used for pairsam parse, migrated from pairtools/pairtools_parse.py
"""

from . import _pairsam_format


def streaming_classify(
    instream,
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
    **kwargs,
):
    """
    TODO: Add handler for pairs parser (regular or complex):
    parser_handler = _parse2.parse_sams_into_pair
    """
    chrom_enum = dict(
        zip(
            [_pairsam_format.UNMAPPED_CHROM] + list(chromosomes),
            range(len(chromosomes) + 1),
        )
    )
    sam_tags = [col for col in add_columns if len(col) == 2 and col.isupper()]
    prev_readID = ""
    sams1 = []
    sams2 = []
    line = ""
    store_seq = "seq" in add_columns

    readID_transform = kwargs.get("readid_transform", None)
    if readID_transform is not None:
        readID_transform = compile(readID_transform, "<string>", "eval")

    instream = iter(instream)
    while line is not None:
        line = next(instream, None)

        readID = line.split("\t", 1)[0] if line else None
        if readID_transform is not None and readID is not None:
            readID = eval(readID_transform)

        if not (line) or ((readID != prev_readID) and prev_readID):

            for (
                algn1,
                algn2,
                all_algns1,
                all_algns2,
                junction_index,
            ) in parse_sams_into_pair(
                sams1,
                sams2,
                min_mapq,
                kwargs["max_inter_align_gap"],
                kwargs["max_fragment_size"],
                sam_tags,
                store_seq,
                kwargs["single_end"],
                coordinate_system,
            ):
                if kwargs["report_alignment_end"] == "5":
                    algn1["pos"] = algn1["pos5"]
                    algn2["pos"] = algn2["pos5"]
                else:
                    algn1["pos"] = algn1["pos3"]
                    algn2["pos"] = algn2["pos3"]


                flip_pair = (not kwargs["no_flip"]) and (
                    not check_pair_order(algn1, algn2, chrom_enum)
                )

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
                    drop_readid,
                    drop_sam,
                    add_junction_index,
                    add_columns,
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

            sams1.clear()
            sams2.clear()

        if line is not None:
            push_sam(line, drop_seq, sams1, sams2)
            prev_readID = readID


def parse_sams_into_pair(
    sams1,
    sams2,
    min_mapq,
    max_inter_align_gap,
    max_fragment_size,
    sam_tags,
    store_seq,
    single_end,
    coordinate_system,
):
    """
    Parse sam entries corresponding to a Hi-C molecule into alignments
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
        sams = sams2  # TODO: Check why it is always the second alignment, and not the first one
        # Generate a sorted, gap-filled list of all alignments
        algns1 = [
            parse_algn(sam.rstrip().split("\t"), min_mapq, sam_tags, store_seq)
            for sam in sams
        ]
        algns1 = sorted(algns1, key=lambda algn: algn["dist_to_5"])
        if max_inter_align_gap is not None:
            _convert_gaps_into_alignments(algns1, max_inter_align_gap)

        algns2 = [empty_alignment()]  # Empty alignment dummy

        if len(algns1) > 1:
            # Look for ligation junction, and report linear alignments after deduplication of complex walks:
            # (Note that coordinate system for single-end reads does not change the behavior)
            return parse_complex_walk(
                algns1, algns2, max_fragment_size, coordinate_system
            )  # TODO: Add offset as param
        else:
            # If no additional information, we assume each molecule is a single ligation with single unconfirmed junction:
            if coordinate_system == "walk":
                return [[algns1[0], flip_alignment(algns2[0]), algns1, algns2, "1u"]]
            else:
                return [[algns1[0], algns2[0], algns1, algns2, "1u"]]

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
        algns1 = [
            parse_algn(sam.rstrip().split("\t"), min_mapq, sam_tags, store_seq)
            for sam in sams1
        ]
        algns2 = [
            parse_algn(sam.rstrip().split("\t"), min_mapq, sam_tags, store_seq)
            for sam in sams2
        ]
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
                algns1, algns2, max_fragment_size, coordinate_system
            )
        else:
            # If no additional information, we assume each molecule is a single ligation with single unconfirmed junction:
            if coordinate_system == "walk":
                return [[algns1[0], flip_alignment(algns2[0]), algns1, algns2, "1u"]]
            else:
                return [[algns1[0], algns2[0], algns1, algns2, "1u"]]


def parse_cigar(cigar):
    matched_bp = 0
    algn_ref_span = 0
    algn_read_span = 0
    read_len = 0
    clip5_ref = 0
    clip3_ref = 0

    if cigar != "*":
        cur_num = 0
        for char in cigar:
            charval = ord(char)
            if charval >= 48 and charval <= 57:
                cur_num = cur_num * 10 + (charval - 48)
            else:
                if char == "M":
                    matched_bp += cur_num
                    algn_ref_span += cur_num
                    algn_read_span += cur_num
                    read_len += cur_num
                elif char == "I":
                    algn_read_span += cur_num
                    read_len += cur_num
                elif char == "D":
                    algn_ref_span += cur_num
                elif char == "S" or char == "H":
                    read_len += cur_num
                    if matched_bp == 0:
                        clip5_ref = cur_num
                    else:
                        clip3_ref = cur_num

                cur_num = 0

    return {
        "clip5_ref": clip5_ref,
        "clip3_ref": clip3_ref,
        "cigar": cigar,
        "algn_ref_span": algn_ref_span,
        "algn_read_span": algn_read_span,
        "read_len": read_len,
        "matched_bp": matched_bp,
    }


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


def parse_algn(samcols, min_mapq, sam_tags=None, store_seq=False):
    is_mapped = (int(samcols[1]) & 0x04) == 0
    mapq = int(samcols[4])
    is_unique = mapq >= min_mapq
    is_linear = not any([col.startswith("SA:Z:") for col in samcols[11:]])

    cigar = parse_cigar(samcols[5])

    if is_mapped:
        if (int(samcols[1]) & 0x10) == 0:
            strand = "+"
            dist_to_5 = cigar["clip5_ref"]
            dist_to_3 = cigar["clip3_ref"]
        else:
            strand = "-"
            dist_to_5 = cigar["clip3_ref"]
            dist_to_3 = cigar["clip5_ref"]

        if is_unique:
            chrom = samcols[2]
            if strand == "+":
                pos5 = int(samcols[3])
                pos3 = int(samcols[3]) + cigar["algn_ref_span"] - 1
            else:
                pos5 = int(samcols[3]) + cigar["algn_ref_span"] - 1
                pos3 = int(samcols[3])
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

    algn["pos"] = algn["pos5"]

    if sam_tags:
        for tag in sam_tags:
            algn[tag] = ""

        for col in samcols[11:]:
            for tag in sam_tags:
                if col.startswith(tag + ":"):
                    algn[tag] = col[5:]
                    continue

    if store_seq:
        algn["seq"] = samcols[9]

    return algn


def parse_complex_walk(
    algns1, algns2, max_fragment_size, coordinate_system, allowed_offset=3
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
    :param coordinate_system: 'walk', 'read' or 'pair'.
        'walk' Alignments are oriented so that 5'-end of each alignment is closer to 5'-end of the walk
        'read' Alignments are oriented so that 5'-end of each alignment is closer to 5'-end of the read
        'pair' Alignments are oriented so that 5'-end of the first alignment is closer to 5'-end of the walk,
            and 5'-end of the second alignment is closer to the 3'-end of the walk
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

    AVAILABLE_COORD_SYSTEMS = ["read", "walk", "pair"]
    assert coordinate_system in AVAILABLE_COORD_SYSTEMS, (
        f"Coordinate system {coordinate_system} is not implemented"
        f'Available choices are: {AVAILABLE_COORD_SYSTEMS.join(", ")}'
    )
    n_algns1 = len(algns1)
    n_algns2 = len(algns2)

    # Iterative search of overlap, let's initialize some useful variables:
    current_forward_junction = current_reverse_junction = 1  # p. 1, initialization
    remaining_forward_junctions = (
        n_algns1 - 1
    )  # Number of possible junctions remaining on forward read
    remaining_reverse_junctions = (
        n_algns2 - 1
    )  # Number of possible junctions remaining on reverse read
    checked_reverse_junctions = (
        0  # Number of checked junctions on reverse read (from the end of read)
    )
    is_overlap = False

    final_contacts = []

    # If both sides have more than 2 alignments, then check if there are overlapping forward and reverse alignments pairs:
    if (n_algns1 >= 2) and (n_algns2 >= 2):
        # Loop through all alignment pairs and check for overlaps:
        while (remaining_forward_junctions > checked_reverse_junctions) and (
            remaining_reverse_junctions > 0
        ):

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
                while is_overlap and (
                    checked_reverse_temp > 0
                ):  # loop over all forward downstream and reverse upstream junctions
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
                if (
                    is_overlap
                ):  # all the checks have passed, no need to check for another hit:
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
            if (
                n_algns1 >= 2
            ):  # Multiple alignments on forward read and single alignment on reverse
                push_pair(
                    algns1[-2],
                    algns2[-1],
                    final_contacts,
                    algns1,
                    algns2,
                    junction_index=f"{len(algns1)-1}f",
                    adjust_reverse_3=algns1[-1][
                        "pos5"
                    ],  # Modify pos3 to correspond to the overlap end at the opposite side
                    flip_reverse=True if coordinate_system == "walk" else False,
                )
                last_reported_alignment_forward = 2
            if (
                n_algns2 >= 2
            ):  # Single alignment on forward read and multiple alignments on reverse
                if coordinate_system == "read":
                    push_pair(
                        algns1[-1],
                        algns2[-2],
                        final_contacts,
                        algns1,
                        algns2,
                        junction_index=f"{len(algns1)}r",
                        adjust_forward_3=algns2[-1][
                            "pos5"
                        ],  # Modify pos3 to correspond to the overlap end at the opposite side
                        flip_forward=True if coordinate_system == "read" else False,
                        flip_reverse=True if coordinate_system == "walk" else False,
                    )
                last_reported_alignment_reverse = 2
            # If n_algns1==n_algns2==1 and alignments overlap, then we don't need to check,
            # it's a non-ligated DNA fragment that we don't report. TODO: rethink this decision?
        # If end alignments do not overlap, then there is no evidence of ligation junction for the pair.
        # Report regular pair:
        else:
            push_pair(
                algns1[-1],
                algns2[-1],
                final_contacts,
                algns1,
                algns2,
                junction_index=f"{len(algns1)}u",
                flip_reverse=True if coordinate_system == "walk" else False,
            )

    # If we have an overlap of junctions:
    else:
        last_reported_alignment_forward = (
            last_reported_alignment_reverse
        ) = current_reverse_junction

    # Report all the sequential alignments:

    # Report all the sequential chimeric pairs in the forward read up to overlap:
    for i in range(0, n_algns1 - last_reported_alignment_forward):
        push_pair(
            algns1[i],
            algns1[i + 1],
            final_contacts,
            algns1,
            algns2,
            junction_index=f"{i + 1}f",
            flip_reverse=True if coordinate_system == "pair" else False,
        )

    # Report the overlap
    for i_overlapping in range(current_reverse_junction - 1):
        idx_forward = n_algns1 - current_reverse_junction + i_overlapping
        idx_reverse = n_algns2 - 1 - i_overlapping
        push_pair(
            algns1[idx_forward],
            algns1[idx_forward + 1],
            final_contacts,
            algns1,
            algns2,
            junction_index=f"{idx_forward + 1}b",
            adjust_reverse_3=algns2[idx_reverse - 1]["pos5"],
            flip_reverse=True if coordinate_system == "pair" else False,
        )

    # Report all the sequential chimeric pairs in the reverse read, but not the overlap:
    for i in range(
        0, min(current_reverse_junction, n_algns2 - last_reported_alignment_reverse)
    ):
        # Two principal cases to determine the junction index for alignments on the reverse read:
        if current_reverse_junction > 1:
            junction_index = (
                n_algns1
                + min(
                    current_reverse_junction, n_algns2 - last_reported_alignment_reverse
                )
                - i
                - 1
            )
        else:
            junction_index = (
                n_algns1
                + min(
                    current_reverse_junction, n_algns2 - last_reported_alignment_reverse
                )
                - i
            )
        if coordinate_system == "pair":
            push_pair(
                algns2[i + 1],
                algns2[i],
                final_contacts,
                algns1,
                algns2,
                junction_index=f"{junction_index}r",
                flip_forward=True,
            )
        elif coordinate_system == "walk":
            push_pair(
                algns2[i + 1],
                algns2[i],
                final_contacts,
                algns1,
                algns2,
                junction_index=f"{junction_index}r",
                flip_forward=True,
                flip_reverse=True,
            )
        else:  # 'read'
            push_pair(
                algns2[i + 1],
                algns2[i],
                final_contacts,
                algns1,
                algns2,
                junction_index=f"{junction_index}r",
            )

    # Sort the pairs according to the order of appearance in the reads.
    # Take the junction index (last element in each entry from its end),
    # and put forward reads first, then the reverse reads:
    final_contacts.sort(key=lambda x: int(x[-1][:-1]))
    return final_contacts


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


def push_pair(
    hic_algn1,
    hic_algn2,
    final_contacts,
    algns1,
    algns2,
    junction_index="1u",
    adjust_forward_3=None,
    adjust_reverse_3=None,
    flip_forward=False,
    flip_reverse=False,
):
    """
    Push a pair of alignments into final list of contacts.
    :param hic_algn1: First alignment in a pair
    :param hic_algn2: Second alignment in a pair
    :param final_contacts: List that will be updated
    :param algns1: All forward read alignments
    :param algns2: All reverse read alignments
    :param junction_index: Index of the junction
    :param adjust_forward_3: Replace 3'-end of the alignment 1 with this position
    :param adjust_forward_3: Replace 3'-end of the alignment 2 with this position
    :return: 0 if successful
    """

    hic_algn1, hic_algn2 = (
        dict(hic_algn1),
        dict(hic_algn2),
    )  # overwrite the variables with copies of dictionaries

    if adjust_forward_3 is not None:  # Adjust forward 3'-end
        hic_algn1["pos3"] = adjust_forward_3
    if adjust_reverse_3 is not None:  # Adjust reverse 3'-end
        hic_algn2["pos3"] = adjust_reverse_3

    hic_algn1["type"] = (
        "N"
        if not hic_algn1["is_mapped"]
        else ("M" if not hic_algn1["is_unique"] else "U")
    )
    hic_algn2["type"] = (
        "N"
        if not hic_algn2["is_mapped"]
        else ("M" if not hic_algn2["is_unique"] else "U")
    )

    if flip_forward:
        hic_algn1 = flip_alignment(hic_algn1)
    if flip_reverse:
        hic_algn2 = flip_alignment(hic_algn2)

    final_contacts.append([hic_algn1, hic_algn2, algns1, algns2, junction_index])

    return 0


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


def _mask_alignment(algn):
    """
    Reset the coordinates of an alignment.
    """
    algn["chrom"] = _pairsam_format.UNMAPPED_CHROM
    algn["pos5"] = _pairsam_format.UNMAPPED_POS
    algn["pos3"] = _pairsam_format.UNMAPPED_POS
    algn["pos"] = _pairsam_format.UNMAPPED_POS
    algn["strand"] = _pairsam_format.UNMAPPED_STRAND

    return algn


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


def push_sam(line, drop_seq, sams1, sams2):
    """
    """

    sam = line.rstrip()
    if drop_seq:
        split_sam = sam.split("\t")
        split_sam[9] = "*"
        split_sam[10] = "*"
        sam = "\t".join(split_sam)

        flag = split_sam[1]
        flag = int(flag)
    else:
        _, flag, _ = sam.split("\t", 2)
        flag = int(flag)

    if (flag & 0x40) != 0:
        sams1.append(sam)
    else:
        sams2.append(sam)
    return


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
            cols.append(
                _pairsam_format.INTER_SAM_SEP.join(
                    [
                        (
                            sam.replace("\t", _pairsam_format.SAM_SEP)
                            + _pairsam_format.SAM_SEP
                            + "Yt:Z:"
                            + algn1["type"]
                            + algn2["type"]
                        )
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
