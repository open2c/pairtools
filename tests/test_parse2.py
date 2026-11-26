# -*- coding: utf-8 -*-
import os
import sys

import pytest

import subprocess

testdir = os.path.dirname(os.path.realpath(__file__))


def test_mock_pysam_parse2_read():
    mock_sam_path = os.path.join(testdir, "data", "mock.parse2.sam")
    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "parse2",
                "-c",
                mock_chroms_path,
                "--add-pair-index",
                "--flip",
                "--report-position",
                "junction",
                "--report-orientation",
                "pair",
                mock_sam_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    # check if the header got transferred correctly
    sam_header = [l.strip() for l in open(mock_sam_path, "r") if l.startswith("@")]
    pairsam_header = [l.strip() for l in result.split("\n") if l.startswith("#")]
    for l in sam_header:
        assert any([l in l2 for l2 in pairsam_header])

    # check that the pairs got assigned properly
    id_counter = 0
    prev_id = ""
    for l in result.split("\n"):
        if l.startswith("#") or not l:
            continue

        if prev_id == l.split("\t")[0]:
            id_counter += 1
        else:
            id_counter = 0
        prev_id = l.split("\t")[0]

        assigned_pair = l.split("\t")[1:8] + l.split("\t")[-2:]
        print(l.split("SIMULATED:", 1)[1].split("\031", 1)[0].split("|"), id_counter)
        simulated_pair = (
            l.split("SIMULATED:", 1)[1]
            .split("\031", 1)[0]
            .split("|")[id_counter]
            .split(",")
        )
        print(assigned_pair)
        print(simulated_pair, prev_id)
        print()

        assert assigned_pair == simulated_pair


def test_mock_pysam_parse2_pair():
    mock_sam_path = os.path.join(testdir, "data", "mock.parse-all.sam")
    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "parse2",
                "-c",
                mock_chroms_path,
                "--add-pair-index",
                "--flip",
                "--report-position",
                "outer",
                "--report-orientation",
                "pair",
                mock_sam_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    # check if the header got transferred correctly
    sam_header = [l.strip() for l in open(mock_sam_path, "r") if l.startswith("@")]
    pairsam_header = [l.strip() for l in result.split("\n") if l.startswith("#")]
    for l in sam_header:
        assert any([l in l2 for l2 in pairsam_header])

    # check that the pairs got assigned properly
    id_counter = 0
    prev_id = ""
    for l in result.split("\n"):
        if l.startswith("#") or not l:
            continue

        if prev_id == l.split("\t")[0]:
            id_counter += 1
        else:
            id_counter = 0
        prev_id = l.split("\t")[0]

        assigned_pair = l.split("\t")[1:8] + l.split("\t")[-2:]
        simulated_pair = (
            l.split("SIMULATED:", 1)[1]
            .split("\031", 1)[0]
            .split("|")[id_counter]
            .split(",")
        )
        print(assigned_pair)
        print(simulated_pair, prev_id)
        print()

        assert assigned_pair == simulated_pair


def test_mock_pysam_parse2_single_end():
    
    """Testing single-end mode for parse2, no-flip mode.
    --report-position is outer (parse2 default)
    --report-orientation is pair (parse2 default)
    """

    mock_sam_path = os.path.join(testdir, "data", "mock.parse2-single-end.sam")
    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "parse2",
                "-c",
                mock_chroms_path,
                "--single-end",
                "--add-pair-index",
                "--no-flip",
                "--report-position",
                "outer",
                "--report-orientation",
                "pair",
                mock_sam_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    # check if the header got transferred correctly
    sam_header = [l.strip() for l in open(mock_sam_path, "r") if l.startswith("@")]
    pairsam_header = [l.strip() for l in result.split("\n") if l.startswith("#")]
    for l in sam_header:
        assert any([l in l2 for l2 in pairsam_header])

    # check that the pairs got assigned properly
    id_counter = 0
    prev_id = ""
    for l in result.split("\n"):
        if l.startswith("#") or not l:
            continue

        if prev_id == l.split("\t")[0]:
            id_counter += 1
        else:
            id_counter = 0
        prev_id = l.split("\t")[0]

        assigned_pair = l.split("\t")[1:8] + l.split("\t")[-2:]
        print(l.split("SIMULATED:", 1)[1].split("\031", 1)[0].split("|"), id_counter)
        simulated_pair = (
            l.split("SIMULATED:", 1)[1]
            .split("\031", 1)[0]
            .split("|")[id_counter]
            .split(",")
        )
        print(assigned_pair)
        print(simulated_pair, prev_id)
        print()

        assert assigned_pair == simulated_pair


def test_mock_pysam_parse2_single_end_expand():
    
    """Testing single-end mode for parse2, no-flip mode, with --expand.
    --report-position is outer (parse2 default)
    --report-orientation is pair (parse2 default)
    """

    mock_sam_path = os.path.join(testdir, "data", "mock.parse2-single-end.expand.sam")
    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "parse2",
                "-c",
                mock_chroms_path,
                "--single-end",
                "--expand",
                "--add-pair-index",
                "--no-flip",
                "--report-position",
                "outer",
                "--report-orientation",
                "pair",
                mock_sam_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    # check if the header got transferred correctly
    sam_header = [l.strip() for l in open(mock_sam_path, "r") if l.startswith("@")]
    pairsam_header = [l.strip() for l in result.split("\n") if l.startswith("#")]
    for l in sam_header:
        assert any([l in l2 for l2 in pairsam_header])

    # check that the pairs got assigned properly
    id_counter = 0
    prev_id = ""
    for l in result.split("\n"):
        if l.startswith("#") or not l:
            continue

        if prev_id == l.split("\t")[0]:
            id_counter += 1
        else:
            id_counter = 0
        prev_id = l.split("\t")[0]

        assigned_pair = l.split("\t")[1:8] + l.split("\t")[-2:]
        print(l.split("SIMULATED:", 1)[1].split("\031", 1)[0].split("|"), id_counter)
        simulated_pair = (
            l.split("SIMULATED:", 1)[1]
            .split("\031", 1)[0]
            .split("|")[id_counter]
            .split(",")
        )
        print(assigned_pair)
        print(simulated_pair, prev_id)
        print()

        assert assigned_pair == simulated_pair



def _write_sam_and_bam(tmp_path, sam_text):
    """Write SAM text to disk and convert it to a BAM file via pysam."""
    import pysam

    sam_path = tmp_path / "input.sam"
    bam_path = tmp_path / "input.bam"
    sam_path.write_text(sam_text)
    with pysam.AlignmentFile(str(sam_path), "r") as fin:
        with pysam.AlignmentFile(str(bam_path), "wb", header=fin.header) as fout:
            for rec in fin:
                fout.write(rec)
    return str(sam_path), str(bam_path)


def _parse_pairs_output(stdout_text):
    """Return (columns_list, list_of_row_dicts) from pairs/pairsam output."""
    columns = None
    rows = []
    for line in stdout_text.splitlines():
        if line.startswith("#columns:"):
            columns = line.split(":", 1)[1].strip().split()
        elif line.startswith("#") or not line:
            continue
        else:
            fields = line.split("\t")
            assert columns is not None
            assert len(fields) == len(columns)
            rows.append(dict(zip(columns, fields)))
    return columns, rows


def test_parse2_add_columns_seq_sam_bam_and_parity(tmp_path):
    """Regression test for #281: parse2 --add-columns seq on SAM and BAM,
    strand orientation, flip swapping, and parity with parse."""
    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")

    seq_fwd_r1 = "ACGTTGCAAACGTTGCAAACGTTGCAAACGTTGCAAACGTTGCAAACGTT"  # 50bp
    seq_fwd_r2 = "GGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGG"  # 50bp
    seq_rev_r1 = "TTAACCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAA"  # 50bp, flag & 16
    seq_rev_r2 = "CCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACC"  # 50bp
    qual50 = "I" * 50

    sam_text = "\n".join(
        [
            "@SQ\tSN:chr1\tLN:10000",
            "@SQ\tSN:chr2\tLN:10000",
            # Pair 1: simple forward-forward on chr1 (no flip needed)
            f"r_simple\t65\tchr1\t10\t60\t50M\tchr1\t200\t0\t{seq_fwd_r1}\t{qual50}\tNM:i:0",
            f"r_simple\t129\tchr1\t200\t60\t50M\tchr1\t10\t0\t{seq_fwd_r2}\t{qual50}\tNM:i:0",
            # Pair 2: R1 on chr2 reverse strand (flag=81), R2 on chr1 forward strand (flag=161)
            # Requires flipping under --flip (chr1 < chr2)
            f"r_flip_rev\t81\tchr2\t200\t60\t50M\tchr1\t10\t0\t{seq_rev_r1}\t{qual50}\tNM:i:0",
            f"r_flip_rev\t161\tchr1\t10\t60\t50M\tchr2\t200\t0\t{seq_rev_r2}\t{qual50}\tNM:i:0",
        ]
    ) + "\n"

    sam_path, bam_path = _write_sam_and_bam(tmp_path, sam_text)

    for input_path in (sam_path, bam_path):
        # 1) parse2 --no-flip (default): sequenced order preserved
        out_noflip = subprocess.check_output(
            [
                sys.executable,
                "-m",
                "pairtools",
                "parse2",
                "-c",
                mock_chroms_path,
                "--no-flip",
                "--drop-sam",
                "--add-columns",
                "seq",
                input_path,
            ]
        ).decode("ascii")
        assert "NotImplemented" not in out_noflip
        cols, rows = _parse_pairs_output(out_noflip)
        assert "seq1" in cols and "seq2" in cols
        by_id = {r["readID"]: r for r in rows}

        assert by_id["r_simple"]["chrom1"] == "chr1"
        assert by_id["r_simple"]["chrom2"] == "chr1"
        assert by_id["r_simple"]["seq1"] == seq_fwd_r1
        assert by_id["r_simple"]["seq2"] == seq_fwd_r2

        # Reverse-strand read preserves exact SAM SEQ verbatim (reference strand)
        assert by_id["r_flip_rev"]["chrom1"] == "chr2"
        assert by_id["r_flip_rev"]["strand1"] == "-"
        assert by_id["r_flip_rev"]["seq1"] == seq_rev_r1
        assert by_id["r_flip_rev"]["chrom2"] == "chr1"
        assert by_id["r_flip_rev"]["strand2"] == "+"
        assert by_id["r_flip_rev"]["seq2"] == seq_rev_r2

        # 2) parse2 --flip: seq1 and seq2 swap along with chrom1/chrom2
        out_flip = subprocess.check_output(
            [
                sys.executable,
                "-m",
                "pairtools",
                "parse2",
                "-c",
                mock_chroms_path,
                "--flip",
                "--drop-sam",
                "--add-columns",
                "seq",
                input_path,
            ]
        ).decode("ascii")
        assert "NotImplemented" not in out_flip
        _, rows_flip = _parse_pairs_output(out_flip)
        by_id_flip = {r["readID"]: r for r in rows_flip}

        assert by_id_flip["r_flip_rev"]["chrom1"] == "chr1"
        assert by_id_flip["r_flip_rev"]["strand1"] == "+"
        assert by_id_flip["r_flip_rev"]["seq1"] == seq_rev_r2
        assert by_id_flip["r_flip_rev"]["chrom2"] == "chr2"
        assert by_id_flip["r_flip_rev"]["strand2"] == "-"
        assert by_id_flip["r_flip_rev"]["seq2"] == seq_rev_r1

        # 3) Parity with pairtools parse --add-columns seq (both flipped)
        out_parse = subprocess.check_output(
            [
                sys.executable,
                "-m",
                "pairtools",
                "parse",
                "-c",
                mock_chroms_path,
                "--flip",
                "--drop-sam",
                "--add-columns",
                "seq",
                input_path,
            ]
        ).decode("ascii")
        assert "NotImplemented" not in out_parse
        _, rows_parse = _parse_pairs_output(out_parse)
        assert rows_flip == rows_parse


def test_parse2_add_columns_seq_missing_and_walks(tmp_path):
    """Missing SEQ ('*' / query_sequence=None) in unmapped reads and
    supplementary walk alignments must emit '*', never 'None' or NotImplementedType."""
    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")

    seq_r1_prim = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"  # 50bp
    seq_r2 = "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG"  # 50bp
    qual50 = "I" * 50

    sam_text = "\n".join(
        [
            "@SQ\tSN:chr1\tLN:10000",
            "@SQ\tSN:chr2\tLN:10000",
            # Unmapped pair where both SEQ fields are '*'
            "r_unmapped_noseq\t77\t*\t0\t0\t*\t*\t0\t0\t*\t*",
            "r_unmapped_noseq\t141\t*\t0\t0\t*\t*\t0\t0\t*\t*",
            # Chimeric walk where primary R1 has SEQ, supplementary R1 has SEQ='*', R2 has SEQ
            f"r_walk_supp_noseq\t65\tchr1\t10\t60\t25M25S\tchr1\t200\t0\t{seq_r1_prim}\t{qual50}\tSA:Z:chr1,5300,-,25M25H,60,0;",
            "r_walk_supp_noseq\t2129\tchr1\t5300\t60\t25M25H\tchr1\t200\t0\t*\t*\tSA:Z:chr1,10,+,25M25S,60,0;",
            f"r_walk_supp_noseq\t129\tchr1\t200\t60\t50M\tchr1\t10\t0\t{seq_r2}\t{qual50}",
            # Missing mate (XX -> empty_alignment on side 1)
            f"r_missing_mate\t129\tchr1\t200\t60\t50M\tchr1\t10\t0\t{seq_r2}\t{qual50}",
        ]
    ) + "\n"

    sam_path, bam_path = _write_sam_and_bam(tmp_path, sam_text)

    for input_path in (sam_path, bam_path):
        out = subprocess.check_output(
            [
                sys.executable,
                "-m",
                "pairtools",
                "parse2",
                "-c",
                mock_chroms_path,
                "--no-flip",
                "--drop-sam",
                "--add-pair-index",
                "--add-columns",
                "seq",
                input_path,
            ]
        ).decode("ascii")
        assert "NotImplemented" not in out
        _, rows = _parse_pairs_output(out)

        unmapped = [r for r in rows if r["readID"] == "r_unmapped_noseq"]
        assert len(unmapped) == 1
        assert unmapped[0]["seq1"] == "*"
        assert unmapped[0]["seq2"] == "*"

        walk_rows = [r for r in rows if r["readID"] == "r_walk_supp_noseq"]
        assert len(walk_rows) == 2
        # Pair 1 (R1): primary R1 (has seq) -> supplementary R1 (SEQ='*')
        assert walk_rows[0]["walk_pair_type"] == "R1"
        assert walk_rows[0]["seq1"] == seq_r1_prim
        assert walk_rows[0]["seq2"] == "*"
        # Pair 2 (R1-2): supplementary R1 (SEQ='*') -> R2 (has seq)
        assert walk_rows[1]["walk_pair_type"] == "R1-2"
        assert walk_rows[1]["seq1"] == "*"
        assert walk_rows[1]["seq2"] == seq_r2

        missing_mate = [r for r in rows if r["readID"] == "r_missing_mate"]
        assert len(missing_mate) == 1
        assert missing_mate[0]["seq1"] not in ("None", "<class 'NotImplementedType'>")
        assert missing_mate[0]["seq2"] not in ("None", "<class 'NotImplementedType'>")


def test_parse2_all_extra_columns_and_drop_seq(tmp_path):
    """Audit all EXTRA_COLUMNS (including seq and mismatches with MD tag) and --drop-seq
    on both SAM and BAM inputs."""
    from pairtools.lib import pairsam_format

    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")

    # 10bp reads; R1 has 1 mismatch at read pos 2 (ref G -> read A), MD:Z:2G7
    seq_r1 = "CCACGTACGT"
    seq_r2 = "TTGGCCAATT"
    qual10 = "IIIIIIIIII"  # Phred 40 ('I' = 73 = 40 + 33)

    sam_text = "\n".join(
        [
            "@SQ\tSN:chr1\tLN:10000",
            "@SQ\tSN:chr2\tLN:10000",
            f"r_all\t65\tchr1\t10\t60\t10M\tchr1\t200\t0\t{seq_r1}\t{qual10}\tNM:i:1\tMD:Z:2G7",
            f"r_all\t129\tchr1\t200\t55\t10M\tchr1\t10\t0\t{seq_r2}\t{qual10}\tNM:i:0\tMD:Z:10",
        ]
    ) + "\n"

    sam_path, bam_path = _write_sam_and_bam(tmp_path, sam_text)
    all_cols = ",".join(pairsam_format.EXTRA_COLUMNS + ["NM"])

    for input_path in (sam_path, bam_path):
        out = subprocess.check_output(
            [
                sys.executable,
                "-m",
                "pairtools",
                "parse2",
                "-c",
                mock_chroms_path,
                "--no-flip",
                "--drop-sam",
                "--add-columns",
                all_cols,
                input_path,
            ]
        ).decode("ascii")
        assert "NotImplemented" not in out
        _, rows = _parse_pairs_output(out)
        assert len(rows) == 1
        row = rows[0]

        assert row["mapq1"] == "60" and row["mapq2"] == "55"
        assert row["pos51"] == "10" and row["pos52"] == "200"
        assert row["pos31"] == "19" and row["pos32"] == "209"
        assert row["cigar1"] == "10M" and row["cigar2"] == "10M"
        assert row["read_len1"] == "10" and row["read_len2"] == "10"
        assert row["matched_bp1"] == "10" and row["matched_bp2"] == "10"
        assert row["algn_ref_span1"] == "10" and row["algn_ref_span2"] == "10"
        assert row["algn_read_span1"] == "10" and row["algn_read_span2"] == "10"
        assert row["dist_to_51"] == "0" and row["dist_to_52"] == "0"
        assert row["dist_to_31"] == "0" and row["dist_to_32"] == "0"
        assert row["seq1"] == seq_r1 and row["seq2"] == seq_r2
        assert row["mismatches1"] == "G:A:40:11:2"
        assert row["mismatches2"] == ""
        assert row["read_side1"] == "1" and row["read_side2"] == "2"
        assert row["algn_idx1"] == "0" and row["algn_idx2"] == "0"
        assert row["same_side_algn_count1"] == "1" and row["same_side_algn_count2"] == "1"
        assert row["NM1"] == "1" and row["NM2"] == "0"

        # --drop-seq with SAM columns kept (exercises query_qualities/query_sequence setters)
        out_drop_seq = subprocess.check_output(
            [
                sys.executable,
                "-m",
                "pairtools",
                "parse2",
                "-c",
                mock_chroms_path,
                "--no-flip",
                "--keep-sam",
                "--drop-seq",
                "--add-columns",
                "seq",
                input_path,
            ]
        ).decode("ascii")
        assert "NotImplemented" not in out_drop_seq
        _, rows_drop = _parse_pairs_output(out_drop_seq)
        assert len(rows_drop) == 1
        # Extracted seq1/seq2 remain intact while sam1/sam2 have SEQ & QUAL stripped to '*'
        assert rows_drop[0]["seq1"] == seq_r1
        assert rows_drop[0]["seq2"] == seq_r2
        assert rows_drop[0]["sam1"].split("\031")[9:11] == ["*", "*"]
        assert rows_drop[0]["sam2"].split("\031")[9:11] == ["*", "*"]
