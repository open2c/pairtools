"""
Set of functions used for pairsam parse, migrated from pairtools/pairtools_parse.py
"""

from . import _pairsam_format


def parse_sams_into_pair(
    sams1,
    sams2,
    min_mapq,
    max_molecule_size,
    max_inter_align_gap,
    walks_policy,
    report_3_alignment_end,
    sam_tags,
    store_seq,
    pysam_backend,
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

    # Check if there is at least one SAM entry per side:
    if (len(sams1) == 0) or (len(sams2) == 0):
        algns1 = [empty_alignment()]
        algns2 = [empty_alignment()]
        algns1[0]["type"] = "X"
        algns2[0]["type"] = "X"
        junction_index = "1u"  # By default, assume each molecule is a single ligation with single unconfirmed junction
        return [[algns1[0], algns2[0], algns1, algns2, junction_index]]

    # Generate a sorted, gap-filled list of all alignments
    algns1 = [
        parse_algn_pysam(sam, min_mapq, report_3_alignment_end, sam_tags, store_seq)
        if pysam_backend
        else parse_algn(
            sam.rstrip().split("\t"),
            min_mapq,
            report_3_alignment_end,
            sam_tags,
            store_seq,
        )
        for sam in sams1
    ]
    algns2 = [
        parse_algn_pysam(sam, min_mapq, report_3_alignment_end, sam_tags, store_seq)
        if pysam_backend
        else parse_algn(
            sam.rstrip().split("\t"),
            min_mapq,
            report_3_alignment_end,
            sam_tags,
            store_seq,
        )
        for sam in sams2
    ]
    algns1 = sorted(algns1, key=lambda algn: algn["dist_to_5"])
    algns2 = sorted(algns2, key=lambda algn: algn["dist_to_5"])

    if max_inter_align_gap is not None:
        _convert_gaps_into_alignments(algns1, max_inter_align_gap)
        _convert_gaps_into_alignments(algns2, max_inter_align_gap)

    # Define the type of alignment on each side.
    # The most important split is between chimeric alignments and linear
    # alignments.

    is_chimeric_1 = len(algns1) > 1
    is_chimeric_2 = len(algns2) > 1

    hic_algn1 = algns1[0]
    hic_algn2 = algns2[0]
    junction_index = "1u"  # By default, assume each molecule is a single ligation with single unconfirmed junction

    # Parse chimeras
    rescued_linear_side = None
    if is_chimeric_1 or is_chimeric_2:

        # Report all the linear alignments in a read pair
        if walks_policy == "all":
            # Report linear alignments after deduplication of complex walks
            return rescue_complex_walk(algns1, algns2, max_molecule_size)

        # Report only two alignments for a read pair
        rescued_linear_side = rescue_walk(algns1, algns2, max_molecule_size)

        # Walk was rescued as a simple walk:
        if rescued_linear_side is not None:
            junction_index = (
                f'{1}{"f" if rescued_linear_side==1 else "r"}'  # TODO: replace
            )
        # Walk is unrescuable:
        else:
            if walks_policy == "mask":
                hic_algn1 = _mask_alignment(dict(hic_algn1))
                hic_algn2 = _mask_alignment(dict(hic_algn2))
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

    return [[hic_algn1, hic_algn2, algns1, algns2, junction_index]]


def parse_cigar(cigar):
    """Parse cigar string."""
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


def parse_cigar_pysam(read):
    """Parse cigar tuples reported as cigartuples of pysam read entry.
    Reports alignment span, clipped nucleotides and more.
    See https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples

    :param read: input pysam read entry
    :return: parsed aligned entry (dictionary)

    """
    matched_bp = 0
    algn_ref_span = 0
    algn_read_span = 0
    read_len = 0
    clip5_ref = 0
    clip3_ref = 0

    cigarstring = read.cigarstring
    cigartuples = read.cigartuples
    if cigartuples is not None:
        for operation, length in cigartuples:
            if operation == 0:  # M, match
                matched_bp += length
                algn_ref_span += length
                algn_read_span += length
                read_len += length
            elif operation == 1:  # I, insertion
                algn_read_span += length
                read_len += length
            elif operation == 2:  # D, deletion
                algn_ref_span += length
            elif (
                operation == 4 or operation == 5
            ):  # S and H, soft clip and hard clip, respectively
                read_len += length
                if matched_bp == 0:
                    clip5_ref = length
                else:
                    clip3_ref = length

    return {
        "clip5_ref": clip5_ref,
        "clip3_ref": clip3_ref,
        "cigar": cigarstring,
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


def parse_algn(
    samcols, min_mapq, report_3_alignment_end=False, sam_tags=None, store_seq=False
):
    """Parse sam alignments."""
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

    algn["pos"] = algn["pos3"] if report_3_alignment_end else algn["pos5"]

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


def parse_algn_pysam(
    sam, min_mapq, report_3_alignment_end=False, sam_tags=None, store_seq=False
):
    """Parse alignments from pysam AlignedSegment entry
    :param sam: input pysam AlignedSegment entry
    :param min_mapq: minimal MAPQ to consider as a proper alignment
    :param report_3_alignment_end: if True, 3'-end of alignment will be reported as position
    :param sam_tags: list of sam tags to store
    :param store_seq: if True, the sequence will be parsed and stored in the output
    :return: parsed aligned entry (dictionary)
    """

    flag = sam.flag
    is_mapped = (flag & 0x04) == 0
    mapq = sam.mapq
    is_unique = sam.is_unique(min_mapq)
    is_linear = sam.is_linear

    cigar = parse_cigar_pysam(sam)

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
                pos5 = (
                    sam.reference_start + 1
                )  # Note that pysam output is zero-based, thus add +1
                pos3 = sam.reference_end + cigar["algn_ref_span"]  # - 1
            else:
                pos5 = sam.reference_start + cigar["algn_ref_span"]  # - 1
                pos3 = (
                    sam.reference_end + 1
                )  # Note that pysam output is zero-based, thus add +1
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


def rescue_complex_walk(algns1, algns2, max_molecule_size, allowed_offset=3):
    """
    Rescue a set of ligations that appear as a complex walk.

    This rescue differs from simple rescue_walk by the step of deduplication.
    If the reads are long enough, the reverse read might read through the forward read's meaningful part.
    If one of the reads contains ligation junction, this might lead to reporting fake contact.
    Thus, the pairs of contacts that overlap are paired-end duplicates and should be reported uniquely.

    Return: list of all the rescued pairs after deduplication with junction index for each pair.

    Example of iterative search (note that it's for the illustration of the algorithm only):

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

    Note that we do not need to perform the shifts on the forward read, because
    biologically overlap can only happen involving both ends of forward and reverse read,
    and shifting one of them is enough.
    """

    n_algns1 = len(algns1)
    n_algns2 = len(algns2)

    # Iterative search of overlap
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

    # If both sides have more than 2 alignments, rescue complex walks
    if (n_algns1 >= 2) and (n_algns2 >= 2):

        # p. 4: if potential overlap can be formed
        while (remaining_forward_junctions > checked_reverse_junctions) and (
            remaining_reverse_junctions > 0
        ):

            # p. 5: check the current pairs of junctions
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

            # p. 5: check the remaining pairs of forward downstream / reverse upstream junctions
            if is_overlap:
                last_idx_forward_temp = current_forward_junction
                last_idx_reverse_temp = current_reverse_junction
                checked_reverse_temp = checked_reverse_junctions
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
                if is_overlap:
                    current_reverse_junction += 1
                    break

            # p. 3: shift the reverse junction pointer by one
            current_reverse_junction += 1
            checked_reverse_junctions += 1
            remaining_reverse_junctions -= 1

        if (
            not is_overlap
        ):  # No overlap found, roll the current_idx_reverse back to the initial value
            current_reverse_junction = 1

    # If no overlapping junctions found, or there are less than 2 chimeras in either forward or reverse read,
    # then current_reverse_junction is 1,
    # check whether the last alignments of forward and reverse reads overlap.
    if current_reverse_junction == 1:
        last_reported_alignment_forward = last_reported_alignment_reverse = 1
        if ends_do_overlap(algns1[-1], algns2[-1], max_molecule_size, allowed_offset):
            # Report the modified last junctions:
            if n_algns1 >= 2:
                # store the type of contact and do not modify original entry:
                hic_algn1 = dict(algns1[-2])
                hic_algn2 = dict(algns2[-1])
                # Modify pos3 of reverse read alignment to correspond to actual observed 5' ends in forward read:
                hic_algn2["pos3"] = algns1[-1]["pos5"]
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
                junction_index = f"{len(algns1)-1}f"
                final_contacts.append(
                    [hic_algn1, hic_algn2, algns1, algns2, junction_index]
                )
                last_reported_alignment_forward = 2
            if n_algns2 >= 2:
                # store the type of contact and do not modify original entry:
                hic_algn1 = dict(algns1[-1])
                hic_algn2 = dict(algns2[-2])
                # Modify pos3 of forward read alignment to correspond to actual observed 5' ends in reverse read:
                hic_algn1["pos3"] = algns2[-1]["pos5"]
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
                junction_index = f"{len(algns1)}r"
                final_contacts.append(
                    [hic_algn1, hic_algn2, algns1, algns2, junction_index]
                )
                last_reported_alignment_reverse = 2
        # End alignments do not overlap. No evidence of ligation junction for the pair, report regular pair:
        else:
            hic_algn1 = dict(
                algns1[-1]
            )  # "dict" trick to store the type of contact and not modify original entry
            hic_algn2 = dict(algns2[-1])
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
            junction_index = f"{len(algns1)}u"
            final_contacts.append(
                [hic_algn1, hic_algn2, algns1, algns2, junction_index]
            )

    # If we have an overlap of junctions:
    else:
        last_reported_alignment_forward = (
            last_reported_alignment_reverse
        ) = current_reverse_junction

    # Report all the sequential alignments
    # Report all the sequential chimeric pairs in the forward read up to overlap:
    for i in range(0, n_algns1 - last_reported_alignment_forward):
        hic_algn1 = dict(algns1[i])
        hic_algn2 = dict(algns1[i + 1])
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
        junction_index = f"{i + 1}f"
        final_contacts.append([hic_algn1, hic_algn2, algns1, algns2, junction_index])

    # Report the overlap
    for i_overlapping in range(current_reverse_junction - 1):
        idx_forward = n_algns1 - current_reverse_junction + i_overlapping
        idx_reverse = n_algns2 - 1 - i_overlapping

        hic_algn1 = dict(algns1[idx_forward])
        hic_algn2 = dict(algns1[idx_forward + 1])
        hic_algn2["pos3"] = algns2[idx_reverse - 1]["pos5"]
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
        junction_index = f"{idx_forward + 1}b"
        final_contacts.append([hic_algn1, hic_algn2, algns1, algns2, junction_index])

    # Report all the sequential chimeric pairs in the reverse read, but not the overlap:
    for i in range(
        0, min(current_reverse_junction, n_algns2 - last_reported_alignment_reverse)
    ):
        hic_algn1 = dict(algns2[i])
        hic_algn2 = dict(algns2[i + 1])
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
        junction_index = f"{n_algns1 +  min(current_reverse_junction, n_algns2 - last_reported_alignment_reverse) - i - (1 if  current_reverse_junction>1 else 0)}r"
        final_contacts.append([hic_algn1, hic_algn2, algns1, algns2, junction_index])

    final_contacts.sort(key=lambda x: int(x[-1][:-1]))
    return final_contacts


### Additional functions for complex walks rescue ###
def ends_do_overlap(algn1, algn2, max_molecule_size=500, allowed_offset=5):
    """
    Two ends of alignments overlap if:
     1) they are from the same chromosome,
     2) map in the opposite directions,
     3) the distance between the outer ends of the two alignments is below the specified max_molecule_size,
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

    do_overlap &= distance_outer_ends <= max_molecule_size + allowed_offset
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
    """Push line into list of sam entries"""

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


def push_pysam(sam, drop_seq, sams1, sams2):
    """Parse pysam AlignedSegment (sam) into pairtools sams entry"""

    flag = sam.flag

    if (flag & 0x40) != 0:
        sams1.append(sam)  # Forward read, or first read in a pair
    else:
        sams2.append(sam)  # Reverse read, or mate pair
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
    pysam_backend,
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
        if pysam_backend:
            for sams in [sams1, sams2]:
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
        else:
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


# TODO: check whether we need this broken function
# def parse_alternative_algns(samcols):
#    alt_algns = []
#    for col in samcols[11:]:
#        if not col.startswith('XA:Z:'):
#            continue
#
#        for SA in col[5:].split(';'):
#            if not SA:
#                continue
#            SAcols = SA.split(',')
#
#            chrom = SAcols[0]
#            strand = '-' if SAcols[1]<0 else '+'
#
#            cigar = parse_cigar(SAcols[2])
#            NM = SAcols[3]
#
#            pos = _pairsam_format.UNMAPPED_POS
#            if strand == '+':
#                pos = int(SAcols[1])
#            else:
#                pos = int(SAcols[1]) + cigar['algn_ref_span']
#
#            alt_algns.append({
#                'chrom': chrom,
#                'pos': pos,
#                'strand': strand,
#                'mapq': mapq, # TODO: Is not defined in this piece of code
#                'is_mapped': True,
#                'is_unique': False,
#                'is_linear': None,
#                'cigar': cigar,
#                'NM': NM,
#                'dist_to_5': cigar['clip5_ref'] if strand == '+' else cigar['clip3_ref'],
#            })
#
#    return supp_algns # TODO: This one seems not to be used in the code...

# def parse_supp(samcols, min_mapq):
#    supp_algns = []
#    for col in samcols[11:]:
#        if not col.startswith('SA:Z:'):
#            continue
#
#        for SA in col[5:].split(';'):
#            if not SA:
#                continue
#            SAcols = SA.split(',')
#            mapq = int(SAcols[4])
#            is_unique = (mapq >= min_mapq)
#
#            chrom = SAcols[0] if is_unique else _pairsam_format.UNMAPPED_CHROM
#            strand = SAcols[2] if is_unique else _pairsam_format.UNMAPPED_STRAND
#
#            cigar = parse_cigar(SAcols[3])
#
#            pos = _pairsam_format.UNMAPPED_POS
#            if is_unique:
#                if strand == '+':
#                    pos = int(SAcols[1])
#                else:
#                    pos = int(SAcols[1]) + cigar['algn_ref_span']
#
#            supp_algns.append({
#                'chrom': chrom,
#                'pos': pos,
#                'strand': strand,
#                'mapq': mapq,
#                'is_mapped': True,
#                'is_unique': is_unique,
#                'is_linear': None,
#                'cigar': cigar,
#                'dist_to_5': cigar['clip5_ref'] if strand == '+' else cigar['clip3_ref'],
#            })
#
#    return supp_algns
