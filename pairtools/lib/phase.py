def get_chrom_phase(chrom, phase_suffixes):
    if chrom.endswith(phase_suffixes[0]):
        return "0", chrom[: -len(phase_suffixes[0])]
    elif chrom.endswith(phase_suffixes[1]):
        return "1", chrom[: -len(phase_suffixes[1])]
    else:
        return "!", chrom


def phase_side_XB(chrom, XB, AS, XS, phase_suffixes):

    phase, chrom_base = get_chrom_phase(chrom, phase_suffixes)

    XBs = [i for i in XB.split(";") if len(i) > 0]
    S1, S2, S3 = AS, XS, -1  # -1 if the second hit was not reported

    if AS > XS:  # Primary hit has higher score than the secondary
        return phase, chrom_base, S1, S2, S3

    elif len(XBs) >= 1:
        if len(XBs) >= 2:
            alt2_chrom, alt2_pos, alt2_CIGAR, alt2_NM, alt2_AS, alt_mapq = XBs[1].split(
                ","
            )
            S3 = int(alt2_AS)
            if int(alt2_AS) == XS == AS:
                return "!", "!", S1, S2, S3

        alt_chrom, alt_pos, alt_CIGAR, alt_NM, alt_AS, alt_mapq = XBs[0].split(",")
        alt_phase, alt_chrom_base = get_chrom_phase(alt_chrom, phase_suffixes)

        alt_is_homologue = (chrom_base == alt_chrom_base) and (
            ((phase == "0") and (alt_phase == "1"))
            or ((phase == "1") and (alt_phase == "0"))
        )

        if alt_is_homologue:
            return ".", chrom_base, S1, S2, S3

    return "!", "!", S1, S2, S3


def phase_side_XA(chrom, XA, AS, XS, NM, phase_suffixes):

    phase, chrom_base = get_chrom_phase(chrom, phase_suffixes)

    XAs = [i for i in XA.split(";") if len(i.strip()) > 0]
    if len(XAs) >= 1:
        alt_chrom, alt_pos, alt_CIGAR, alt_NM = XAs[0].split(",")
        M1, M2, M3 = NM, int(alt_NM), -1
    else:
        M1, M2, M3 = NM, -1, -1  # -1 if the second hit was not reported

    if AS > XS:  # Primary hit has higher score than the secondary
        return phase, chrom_base, M1, M2, M3

    elif len(XAs) >= 1:

        if len(XAs) >= 2:
            alt2_chrom, alt2_pos, alt2_CIGAR, alt2_NM = XAs[1].split(",")
            M3 = int(alt2_NM)
            if int(alt2_NM) == int(alt_NM) == NM:
                return "!", "!", M1, M2, M3

        alt_chrom, alt_pos, alt_CIGAR, alt_NM = XAs[0].split(",")

        alt_phase, alt_chrom_base = get_chrom_phase(alt_chrom, phase_suffixes)

        alt_is_homologue = (chrom_base == alt_chrom_base) and (
            ((phase == "0") and (alt_phase == "1"))
            or ((phase == "1") and (alt_phase == "0"))
        )

        if alt_is_homologue:
            return ".", chrom_base, M1, M2, M3

    return "!", "!", M1, M2, M3
