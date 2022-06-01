from . import pairsam_format
import warnings


def find_rfrag(rfrags, chrom, pos):

    # Return empty if chromosome is unmapped:
    if chrom == pairsam_format.UNMAPPED_CHROM:
        return (
            pairsam_format.UNANNOTATED_RFRAG,
            pairsam_format.UNMAPPED_POS,
            pairsam_format.UNMAPPED_POS,
        )

    try:
        rsites_chrom = rfrags[chrom]
    except ValueError as e:
        warnings.warn(
            f"Chomosome {chrom} does not have annotated restriction fragments, return empty."
        )
        return (
            pairsam_format.UNANNOTATED_RFRAG,
            pairsam_format.UNMAPPED_POS,
            pairsam_format.UNMAPPED_POS,
        )

    idx = min(
        max(0, rsites_chrom.searchsorted(pos, "right") - 1), len(rsites_chrom) - 2
    )
    return idx, rsites_chrom[idx], rsites_chrom[idx + 1]
