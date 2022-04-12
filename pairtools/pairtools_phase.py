import sys
import click
import re, fnmatch

from . import _fileio, _pairsam_format, cli, _headerops, common_io_options

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
    help="phase suffixes.",
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
    " XB will parse 'XB' tag, the input should be generated with: --add-columns XB,AS,XS --min-mapq 0 "
    " XA will parse 'XA', the input should be generated with: --add-columns XA,NM,AS,XS --min-mapq 0"
    " Note that XB tag be added by running bwa with -u tag, present in github version.",
)
@click.option(
    "--report-scores/--no-report-scores",
    is_flag=True,
    default=False,
    help="Report scores of optinal, suboptimal and second suboptimal alignments. "
         "NM with --tag-mode XA and AS with --tag-mode XB ",
)
@common_io_options
def phase(pairs_path, output, phase_suffixes, clean_output, tag_mode, report_scores, **kwargs):
    """Phase pairs mapped to a diploid genome.

    PAIRS_PATH : input .pairs/.pairsam file. If the path ends with .gz or .lz4, the
    input is decompressed by bgzip/lz4c. By default, the input is read from stdin.

    """
    phase_py(pairs_path, output, phase_suffixes, clean_output, tag_mode, report_scores, **kwargs)


if __name__ == "__main__":
    phase()


def phase_py(pairs_path, output, phase_suffixes, clean_output, tag_mode, report_scores, **kwargs):

    instream = (
        _fileio.auto_open(
            pairs_path,
            mode="r",
            nproc=kwargs.get("nproc_in"),
            command=kwargs.get("cmd_in", None),
        )
        if pairs_path
        else sys.stdin
    )
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

    header, body_stream = _headerops.get_header(instream)
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    old_column_names = _headerops.extract_column_names(header)

    idx_phase1 = len(old_column_names)
    idx_phase2 = len(old_column_names) + 1
    if clean_output:
        new_column_names = [
            col for col in old_column_names if col in _pairsam_format.COLUMNS
        ]
        new_column_idxs = [
            i
            for i, col in enumerate(old_column_names)
            if col in _pairsam_format.COLUMNS
        ]
        new_column_idxs += [idx_phase1, idx_phase2]
    else:
        new_column_names = list(old_column_names)

    new_column_names.append("phase1")
    new_column_names.append("phase2")

    if report_scores:
        if tag_mode=="XB":
            new_column_names.append("S1_1")
            new_column_names.append("S1_2")
            new_column_names.append("S2_1")
            new_column_names.append("S2_2")
            new_column_names.append("S3_1")
            new_column_names.append("S3_2")
            if clean_output:
                new_column_idxs += [(idx_phase2 + i + 1) for i in range(6)]
        elif tag_mode=="XA":
            new_column_names.append("M1_1")
            new_column_names.append("M1_2")
            new_column_names.append("M2_1")
            new_column_names.append("M2_2")
            new_column_names.append("M3_1")
            new_column_names.append("M3_2")
            if clean_output:
                new_column_idxs += [(idx_phase2 + i + 1) for i in range(6)]
    header = _headerops._update_header_entry(
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
        cols = line.rstrip().split(_pairsam_format.PAIRSAM_SEP)
        cols.append("!")
        cols.append("!")
        if report_scores:
            for _ in range(6):
                cols.append("!")
        pair_type = cols[_pairsam_format.COL_PTYPE]

        if cols[_pairsam_format.COL_C1] != _pairsam_format.UNMAPPED_CHROM:
            if tag_mode == "XB":
                phase1, chrom_base1, S1_1, S2_1, S3_1 = phase_side_XB(
                    cols[_pairsam_format.COL_C1],
                    cols[COL_XB1],
                    int(cols[COL_AS1]),
                    int(cols[COL_XS1]),
                    phase_suffixes,
                )
            elif tag_mode == "XA":
                phase1, chrom_base1, S1_1, S2_1, S3_1 = phase_side_XA(
                    cols[_pairsam_format.COL_C1],
                    cols[COL_XA1],
                    int(cols[COL_AS1]),
                    int(cols[COL_XS1]),
                    int(cols[COL_NM1]),
                    phase_suffixes,
                )

            if not report_scores:
                cols[idx_phase1] = phase1
            else:
                cols[idx_phase1], cols[idx_phase1+2], cols[idx_phase1+4], cols[idx_phase1+6] \
                    = phase1, str(S1_1), str(S2_1), str(S3_1)
            cols[_pairsam_format.COL_C1] = chrom_base1

            if chrom_base1 == "!":
                cols[_pairsam_format.COL_C1] = _pairsam_format.UNMAPPED_CHROM
                cols[_pairsam_format.COL_P1] = str(_pairsam_format.UNMAPPED_POS)
                cols[_pairsam_format.COL_S1] = _pairsam_format.UNMAPPED_STRAND
                pair_type = "M" + pair_type[1]

        if cols[_pairsam_format.COL_C2] != _pairsam_format.UNMAPPED_CHROM:

            if tag_mode == "XB":
                phase2, chrom_base2, S1_2, S2_2, S3_2 = phase_side_XB(
                    cols[_pairsam_format.COL_C2],
                    cols[COL_XB2],
                    int(cols[COL_AS2]),
                    int(cols[COL_XS2]),
                    phase_suffixes,
                )
            elif tag_mode == "XA":
                phase2, chrom_base2, S1_2, S2_2, S3_2 = phase_side_XA(
                    cols[_pairsam_format.COL_C2],
                    cols[COL_XA2],
                    int(cols[COL_AS2]),
                    int(cols[COL_XS2]),
                    int(cols[COL_NM2]),
                    phase_suffixes,
                )

            if not report_scores:
                cols[idx_phase1] = phase2
            else:
                cols[idx_phase2], cols[idx_phase2+2], cols[idx_phase2+4], cols[idx_phase2+6] \
                    = phase2, str(S1_2), str(S2_2), str(S3_2)
            cols[_pairsam_format.COL_C2] = chrom_base2

            if chrom_base2 == "!":
                cols[_pairsam_format.COL_C2] = _pairsam_format.UNMAPPED_CHROM
                cols[_pairsam_format.COL_P2] = str(_pairsam_format.UNMAPPED_POS)
                cols[_pairsam_format.COL_S2] = _pairsam_format.UNMAPPED_STRAND
                pair_type = pair_type[0] + "M"

        cols[_pairsam_format.COL_PTYPE] = pair_type

        if clean_output:
            cols = [cols[i] for i in new_column_idxs]

        outstream.write(_pairsam_format.PAIRSAM_SEP.join(cols))
        outstream.write("\n")

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()


def get_chrom_phase(chrom, phase_suffixes):
    if chrom.endswith(phase_suffixes[0]):
        return "0", chrom[: -len(phase_suffixes[0])]
    elif chrom.endswith(phase_suffixes[1]):
        return "1", chrom[: -len(phase_suffixes[1])]
    else:
        return "!", chrom


def phase_side_XB(chrom, XB, AS, XS, phase_suffixes):

    phase, chrom_base = get_chrom_phase(chrom, phase_suffixes)

    XBs = [i for i in XB.split(';') if len(i) > 0]
    S1, S2, S3 = AS, XS, -1 # -1 if the second hit was not reported

    if AS > XS: # Primary hit has higher score than the secondary
        return phase, chrom_base, S1, S2, S3

    elif len(XBs) >= 1:
        if len(XBs) >= 2:
            alt2_chrom, alt2_pos, alt2_CIGAR, alt2_NM, alt2_AS, alt_mapq = XBs[1].split(',')
            S3 = int(alt2_AS)
            if int(alt2_AS) == XS == AS:
                return '!', '!', S1, S2, S3

        alt_chrom, alt_pos, alt_CIGAR, alt_NM, alt_AS, alt_mapq = XBs[0].split(',')
        alt_phase, alt_chrom_base = get_chrom_phase(alt_chrom, phase_suffixes)

        alt_is_homologue = (
                (chrom_base == alt_chrom_base)
                and
                (
                        ((phase == '0') and (alt_phase == '1'))
                        or
                        ((phase == '1') and (alt_phase == '0'))
                )
        )

        if alt_is_homologue:
            return '.', chrom_base, S1, S2, S3

    return '!', '!', S1, S2, S3


def phase_side_XA(chrom, XA, AS, XS, NM, phase_suffixes):

    phase, chrom_base = get_chrom_phase(chrom, phase_suffixes)

    XAs = [i for i in XA.split(";") if len(i.strip()) > 0]
    if len(XAs) >= 1:
        alt_chrom, alt_pos, alt_CIGAR, alt_NM = XAs[0].split(",")
        M1, M2, M3 = NM, int(alt_NM), -1
    else:
        M1, M2, M3 = NM, -1, -1 # -1 if the second hit was not reported

    if (AS > XS):  # Primary hit has higher score than the secondary
        return phase, chrom_base, M1, M2, M3

    elif len(XAs) >= 1:

        if len(XAs) >= 2:
            alt2_chrom, alt2_pos, alt2_CIGAR, alt2_NM = XAs[1].split(",")
            M3 = int(alt2_NM)
            if int(alt2_NM) == NM:
                return "!", "!", M1, M2, M3

        alt_chrom, alt_pos, alt_CIGAR, alt_NM = XAs[0].split(",")

        # if NM < int(alt_NM):
        #     return phase, chrom_base

        alt_phase, alt_chrom_base = get_chrom_phase(alt_chrom, phase_suffixes)

        alt_is_homologue = (chrom_base == alt_chrom_base) and (
            ((phase == "0") and (alt_phase == "1"))
            or ((phase == "1") and (alt_phase == "0"))
        )

        if alt_is_homologue:
            return ".", chrom_base, M1, M2, M3

    return "!", "!", M1, M2, M3
