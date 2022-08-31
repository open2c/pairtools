# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import pytest
import tempfile

testdir = os.path.dirname(os.path.realpath(__file__))
mock_pairsam_path = os.path.join(testdir, "data", "mock.pairsam")

tmpdir = tempfile.TemporaryDirectory()
tmpdir_name = tmpdir.name
pairs_path = os.path.join(tmpdir_name, "out.pairs")
sam_path = os.path.join(tmpdir_name, "out.sam")


@pytest.fixture
def setup_split():
    try:
        subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "split",
                mock_pairsam_path,
                "--output-pairs",
                pairs_path,
                "--output-sam",
                sam_path,
            ],
        )
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e


def test_split(setup_split):

    pairsam_lines = [l.strip() for l in open(mock_pairsam_path, "r") if l.strip()]
    pairs_lines = [l.strip() for l in open(pairs_path, "r") if l.strip()]
    sam_lines = [l.strip() for l in open(sam_path, "r") if l.strip()]

    # check that all entries survived splitting:
    n_pairsam = len([l for l in pairsam_lines if not l.startswith("#")])
    n_pairs = len([l for l in pairs_lines if not l.startswith("#")])
    n_sam = len([l for l in sam_lines if not l.startswith("@")]) // 2

    assert n_pairsam == n_pairs
    assert n_pairsam == n_sam

    # check that the header survived splitting:
    pairsam_header = [
        l.strip()
        for l in open(mock_pairsam_path, "r")
        if l.strip() and l.startswith("#")
    ]
    pairs_header = [
        l.strip() for l in open(pairs_path, "r") if l.strip() and l.startswith("#")
    ]
    sam_header = [
        l.strip() for l in open(sam_path, "r") if l.strip() and l.startswith("@")
    ]
    assert all(
        any(l in l2 for l2 in pairsam_header)
        for l in sam_header
        if not l.startswith("@PG")
    )
    assert all(
        l in pairsam_header
        for l in pairs_header
        if (not (l.startswith("#columns") or l.startswith("#samheader")))
    )
    columns_pairsam = [l for l in pairsam_header if l.startswith("#columns")][
        0
    ].split()[1:]
    columns_pairs = [l for l in pairs_header if l.startswith("#columns")][0].split()[1:]
    assert (
        ("sam1" in columns_pairsam)
        and ("sam2" in columns_pairsam)
        and ("sam1" not in columns_pairs)
        and ("sam2" not in columns_pairs)
    )
    assert [c for c in columns_pairsam if c != "sam1" and c != "sam2"] == columns_pairs

    tmpdir.cleanup()
