# -*- coding: utf-8 -*-
import os
import sys

import pytest

import subprocess

testdir = os.path.dirname(os.path.realpath(__file__))


def test_python_version():
    assert sys.version_info[0] == 3, "Use Python 3!"


def test_mock_pysam():
    """Parse non-chimeric alignments with walks-policy mask with pysam backend."""
    mock_sam_path = os.path.join(testdir, "data", "mock.sam")
    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "parse",
                "--walks-policy",
                "mask",
                "-c",
                mock_chroms_path,
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
    for l in result.split("\n"):
        if l.startswith("#") or not l:
            continue

        print(l)
        assigned_pair = l.split("\t")[1:8]
        simulated_pair = l.split("CT:Z:SIMULATED:", 1)[1].split("\031", 1)[0].split(",")
        print(assigned_pair)
        print(simulated_pair)
        print()

        assert assigned_pair == simulated_pair


def test_mock_pysam_parse_all():
    """Parse all alignment in each read with walks-policy all and pysam backend."""
    mock_sam_path = os.path.join(testdir, "data", "mock.parse-all.sam")
    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "parse",
                "--walks-policy",
                "all",
                "-c",
                mock_chroms_path,
                "--add-pair-index",
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
            l.split("CT:Z:SIMULATED:", 1)[1]
            .split("\031", 1)[0]
            .split("|")[id_counter]
            .split(",")
        )
        print(assigned_pair)
        print(simulated_pair, prev_id)
        print()

        assert assigned_pair == simulated_pair



def test_parse_add_columns_seq_sam_and_bam(tmp_path):
    """Regression test for #281 on pairtools parse and AlignmentFilePairtoolized."""
    import pysam
    from pairtools.lib.parse_pysam import AlignmentFilePairtoolized

    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")
    seq1 = "ACGTTGCAAACGTTGCAAACGTTGCAAACGTTGCAAACGTTGCAAACGTT"
    seq2 = "GGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGG"
    qual50 = "I" * 50

    sam_text = "\n".join(
        [
            "@SQ\tSN:chr1\tLN:1000",
            "@SQ\tSN:chr2\tLN:1000",
            f"r1\t65\tchr1\t10\t60\t50M\tchr1\t200\t0\t{seq1}\t{qual50}",
            f"r1\t129\tchr1\t200\t60\t50M\tchr1\t10\t0\t{seq2}\t{qual50}",
            "r2\t77\t*\t0\t0\t*\t*\t0\t0\t*\t*",
            "r2\t141\t*\t0\t0\t*\t*\t0\t0\t*\t*",
        ]
    ) + "\n"

    sam_path = tmp_path / "parse_seq.sam"
    bam_path = tmp_path / "parse_seq.bam"
    sam_path.write_text(sam_text)
    with pysam.AlignmentFile(str(sam_path), "r") as fin:
        with pysam.AlignmentFile(str(bam_path), "wb", header=fin.header) as fout:
            for rec in fin:
                fout.write(rec)

    # Verify AlignmentFilePairtoolized cached properties directly
    for path in (str(sam_path), str(bam_path)):
        with AlignmentFilePairtoolized(path, "r") as af:
            segs = list(af)
        assert segs[0].query_sequence == seq1
        assert segs[0].seq == seq1
        assert segs[0].query_alignment_sequence == seq1
        assert len(segs[0].query_qualities) == 50
        assert segs[2].query_sequence is None

        out = subprocess.check_output(
            [
                sys.executable,
                "-m",
                "pairtools",
                "parse",
                "-c",
                mock_chroms_path,
                "--drop-sam",
                "--add-columns",
                "seq",
                path,
            ]
        ).decode("ascii")
        assert "NotImplemented" not in out
        body = [l.split("\t") for l in out.splitlines() if l and not l.startswith("#")]
        by_id = {row[0]: row for row in body}
        assert by_id["r1"][-2:] == [seq1, seq2]
        assert by_id["r2"][-2:] == ["*", "*"]
