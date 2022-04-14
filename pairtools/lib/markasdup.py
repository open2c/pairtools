from . import pairsam_format


def mark_split_pair_as_dup(cols):
    # if the original columns ended with a new line, the marked columns
    # should as well.
    original_has_newline = cols[-1].endswith("\n")

    cols[pairsam_format.COL_PTYPE] = "DD"

    if (len(cols) > pairsam_format.COL_SAM1) and (len(cols) > pairsam_format.COL_SAM2):
        for i in (pairsam_format.COL_SAM1, pairsam_format.COL_SAM2):

            # split each sam column into sam entries, tag and assemble back
            cols[i] = pairsam_format.INTER_SAM_SEP.join(
                [
                    mark_sam_as_dup(sam)
                    for sam in cols[i].split(pairsam_format.INTER_SAM_SEP)
                ]
            )

    if original_has_newline and not cols[-1].endswith("\n"):
        cols[-1] = cols[-1] + "\n"
    return cols


def mark_sam_as_dup(sam):
    """Tag the binary flag and the optional pair type field of a sam entry
    as a PCR duplicate."""
    samcols = sam.split(pairsam_format.SAM_SEP)

    if len(samcols) == 1:
        return sam

    samcols[1] = str(int(samcols[1]) | 1024)

    for j in range(11, len(samcols)):
        if samcols[j].startswith("Yt:Z:"):
            samcols[j] = "Yt:Z:DD"
    return pairsam_format.SAM_SEP.join(samcols)
