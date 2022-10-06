import sys
import click
import re, fnmatch

from ..lib import fileio, pairsam_format, headerops
from . import cli, common_io_options

from ..lib.phase import phase_side_XB, phase_side_XA

UTIL_NAME = "pairtools_phase"


@cli.command()
@click.argument("pairs_path", type=str, required=False)
@click.option(
    "-o",
    "--output",
    type=str,
    default="",
    help="output file."
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4c-compressed."
    " By default, the output is printed into stdout.",
)
@click.option(
    "--phase-suffixes",
    nargs=2,
    # type=click.Tuple([str, str]),
    help="Phase suffixes (of the chrom names), always a pair.",
)
@click.option(
    "--clean-output",
    is_flag=True,
    help="Drop all columns besides the standard ones and phase1/2",
)
@click.option(
    "--tag-mode",
    type=click.Choice(["XB", "XA"]),
    default="XB",
    help="Specifies the mode of bwa reporting."
    " XA will parse 'XA', the input should be generated with: --add-columns XA,NM,AS,XS --min-mapq 0"
    " XB will parse 'XB' tag, the input should be generated with: --add-columns XB,AS,XS --min-mapq 0 "
    " Note that XB tag is added by running bwa with -u tag, present in github version. "
    " Both modes report similar results: XB reports 0.002% contacts more for phased data, "
    " while XA can report ~1-2% more unphased contacts because its definition multiple mappers is more premissive. ",
)
@click.option(
    "--report-scores/--no-report-scores",
    is_flag=True,
    default=False,
    help="Report scores of optional, suboptimal and second suboptimal alignments. "
    "NM (edit distance) with --tag-mode XA and AS (alfn score) with --tag-mode XB ",
)
@common_io_options
def phase(
    pairs_path, output, phase_suffixes, clean_output, tag_mode, report_scores, **kwargs
):
    """Phase pairs mapped to a diploid genome.
    Diploid genome is the genome with two set of the chromosome variants,
    where each chromosome has one of two suffixes (phase-suffixes)
    corresponding to the genome version (phase-suffixes).

    By default, phasing adds two additional columns with phase 0, 1 or "." (unpahsed).

    Phasing is based on detection of chromosome origin of each mapped fragment.
    Three scores are considered: best alignment score (S1),
    suboptimal alignment (S2) and second suboptimal alignment (S3) scores.
    Each fragment can be:
    1) uniquely mapped and phased (S1>S2>S3, first alignment is the best hit),
    2) uniquely mapped but unphased (S1=S2>S3, cannot distinguish between chromosome variants),
    3) multiply mapped (S1=S2=S3) or unmapped.

    PAIRS_PATH : input .pairs/.pairsam file. If the path ends with .gz or .lz4, the
    input is decompressed by bgzip/lz4c. By default, the input is read from stdin.

    """
    phase_py(
        pairs_path,
        output,
        phase_suffixes,
        clean_output,
        tag_mode,
        report_scores,
        **kwargs
    )


if __name__ == "__main__":
    phase()


def phase_py(
    pairs_path, output, phase_suffixes, clean_output, tag_mode, report_scores, **kwargs
):

    instream = (
        fileio.auto_open(
            pairs_path,
            mode="r",
            nproc=kwargs.get("nproc_in"),
            command=kwargs.get("cmd_in", None),
        )
        if pairs_path
        else sys.stdin
    )
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

    header, body_stream = headerops.get_header(instream)
    header = headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    old_column_names = headerops.extract_column_names(header)

    idx_phase1 = len(old_column_names)
    idx_phase2 = len(old_column_names) + 1
    if clean_output:
        new_column_names = [
            col for col in old_column_names if col in pairsam_format.COLUMNS
        ]
        new_column_idxs = [
            i for i, col in enumerate(old_column_names) if col in pairsam_format.COLUMNS
        ]
        new_column_idxs += [idx_phase1, idx_phase2]
    else:
        new_column_names = list(old_column_names)

    new_column_names.append("phase1")
    new_column_names.append("phase2")

    if report_scores:
        if tag_mode == "XB":
            new_column_names.append("S1_1")
            new_column_names.append("S1_2")
            new_column_names.append("S2_1")
            new_column_names.append("S2_2")
            new_column_names.append("S3_1")
            new_column_names.append("S3_2")
            if clean_output:
                new_column_idxs += [(idx_phase2 + i + 1) for i in range(6)]
        elif tag_mode == "XA":
            new_column_names.append("M1_1")
            new_column_names.append("M1_2")
            new_column_names.append("M2_1")
            new_column_names.append("M2_2")
            new_column_names.append("M3_1")
            new_column_names.append("M3_2")
            if clean_output:
                new_column_idxs += [(idx_phase2 + i + 1) for i in range(6)]
    header = headerops._update_header_entry(
        header, "columns", " ".join(new_column_names)
    )

    if tag_mode == "XB":
        if (
            ("XB1" not in old_column_names)
            or ("XB2" not in old_column_names)
            or ("AS1" not in old_column_names)
            or ("AS2" not in old_column_names)
            or ("XS1" not in old_column_names)
            or ("XS2" not in old_column_names)
        ):
            raise ValueError(
                "The input pairs file must be parsed with the flag --add-columns XB,AS,XS --min-mapq 0"
            )

        COL_XB1 = old_column_names.index("XB1")
        COL_XB2 = old_column_names.index("XB2")
        COL_AS1 = old_column_names.index("AS1")
        COL_AS2 = old_column_names.index("AS2")
        COL_XS1 = old_column_names.index("XS1")
        COL_XS2 = old_column_names.index("XS2")

    elif tag_mode == "XA":
        if (
            ("XA1" not in old_column_names)
            or ("XA2" not in old_column_names)
            or ("NM1" not in old_column_names)
            or ("NM2" not in old_column_names)
            or ("AS1" not in old_column_names)
            or ("AS2" not in old_column_names)
            or ("XS1" not in old_column_names)
            or ("XS2" not in old_column_names)
        ):
            raise ValueError(
                "The input pairs file must be parsed with the flag --add-columns XA,NM,AS,XS --min-mapq 0"
            )

        COL_XA1 = old_column_names.index("XA1")
        COL_XA2 = old_column_names.index("XA2")
        COL_NM1 = old_column_names.index("NM1")
        COL_NM2 = old_column_names.index("NM2")
        COL_AS1 = old_column_names.index("AS1")
        COL_AS2 = old_column_names.index("AS2")
        COL_XS1 = old_column_names.index("XS1")
        COL_XS2 = old_column_names.index("XS2")

    outstream.writelines((l + "\n" for l in header))

    for line in body_stream:
        cols = line.rstrip('\n').split(pairsam_format.PAIRSAM_SEP)
        cols.append("!")
        cols.append("!")
        if report_scores:
            for _ in range(6):
                cols.append("-1")
        pair_type = cols[pairsam_format.COL_PTYPE]

        if cols[pairsam_format.COL_C1] != pairsam_format.UNMAPPED_CHROM:
            if tag_mode == "XB":
                phase1, chrom_base1, S1_1, S2_1, S3_1 = phase_side_XB(
                    cols[pairsam_format.COL_C1],
                    cols[COL_XB1],
                    int(cols[COL_AS1]),
                    int(cols[COL_XS1]),
                    phase_suffixes,
                )
            elif tag_mode == "XA":
                phase1, chrom_base1, S1_1, S2_1, S3_1 = phase_side_XA(
                    cols[pairsam_format.COL_C1],
                    cols[COL_XA1],
                    int(cols[COL_AS1]),
                    int(cols[COL_XS1]),
                    int(cols[COL_NM1]),
                    phase_suffixes,
                )

            if not report_scores:
                cols[idx_phase1] = phase1
            else:
                (
                    cols[idx_phase1],
                    cols[idx_phase1 + 2],
                    cols[idx_phase1 + 4],
                    cols[idx_phase1 + 6],
                ) = (phase1, str(S1_1), str(S2_1), str(S3_1))
            cols[pairsam_format.COL_C1] = chrom_base1

            if chrom_base1 == "!":
                cols[pairsam_format.COL_C1] = pairsam_format.UNMAPPED_CHROM
                cols[pairsam_format.COL_P1] = str(pairsam_format.UNMAPPED_POS)
                cols[pairsam_format.COL_S1] = pairsam_format.UNMAPPED_STRAND
                pair_type = "M" + pair_type[1]

        if cols[pairsam_format.COL_C2] != pairsam_format.UNMAPPED_CHROM:

            if tag_mode == "XB":
                phase2, chrom_base2, S1_2, S2_2, S3_2 = phase_side_XB(
                    cols[pairsam_format.COL_C2],
                    cols[COL_XB2],
                    int(cols[COL_AS2]),
                    int(cols[COL_XS2]),
                    phase_suffixes,
                )
            elif tag_mode == "XA":
                phase2, chrom_base2, S1_2, S2_2, S3_2 = phase_side_XA(
                    cols[pairsam_format.COL_C2],
                    cols[COL_XA2],
                    int(cols[COL_AS2]),
                    int(cols[COL_XS2]),
                    int(cols[COL_NM2]),
                    phase_suffixes,
                )
            if not report_scores:
                cols[idx_phase2] = phase2
            else:
                (
                    cols[idx_phase2],
                    cols[idx_phase2 + 2],
                    cols[idx_phase2 + 4],
                    cols[idx_phase2 + 6],
                ) = (phase2, str(S1_2), str(S2_2), str(S3_2))
            cols[pairsam_format.COL_C2] = chrom_base2

            if chrom_base2 == "!":
                cols[pairsam_format.COL_C2] = pairsam_format.UNMAPPED_CHROM
                cols[pairsam_format.COL_P2] = str(pairsam_format.UNMAPPED_POS)
                cols[pairsam_format.COL_S2] = pairsam_format.UNMAPPED_STRAND
                pair_type = pair_type[0] + "M"

        cols[pairsam_format.COL_PTYPE] = pair_type

        if clean_output:
            cols = [cols[i] for i in new_column_idxs]

        outstream.write(pairsam_format.PAIRSAM_SEP.join(cols))
        outstream.write("\n")

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()


if __name__ == "__main__":
    phase()
