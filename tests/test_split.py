# -*- coding: utf-8 -*-
import os
import subprocess
import sys
import tempfile

import pytest

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

    pairsam_lines = [line.strip() for line in open(mock_pairsam_path, "r") if line.strip()]
    pairs_lines = [line.strip() for line in open(pairs_path, "r") if line.strip()]
    sam_lines = [line.strip() for line in open(sam_path, "r") if line.strip()]

    # check that all entries survived splitting:
    n_pairsam = len([line for line in pairsam_lines if not line.startswith("#")])
    n_pairs = len([line for line in pairs_lines if not line.startswith("#")])
    n_sam = len([line for line in sam_lines if not line.startswith("@")]) // 2

    assert n_pairsam == n_pairs
    assert n_pairsam == n_sam

    # check that the header survived splitting:
    pairsam_header = [
        line.strip()
        for line in open(mock_pairsam_path, "r")
        if line.strip() and line.startswith("#")
    ]
    pairs_header = [
        line.strip() for line in open(pairs_path, "r") if line.strip() and line.startswith("#")
    ]
    sam_header = [
        line.strip() for line in open(sam_path, "r") if line.strip() and line.startswith("@")
    ]
    assert all(
        any(line in l2 for l2 in pairsam_header)
        for line in sam_header
        if not line.startswith("@PG")
    )
    assert all(
        line in pairsam_header
        for line in pairs_header
        if (not (line.startswith("#columns") or line.startswith("#samheader")))
    )
    columns_pairsam = [line for line in pairsam_header if line.startswith("#columns")][
        0
    ].split()[1:]
    columns_pairs = [line for line in pairs_header if line.startswith("#columns")][0].split()[1:]
    assert (
        ("sam1" in columns_pairsam)
        and ("sam2" in columns_pairsam)
        and ("sam1" not in columns_pairs)
        and ("sam2" not in columns_pairs)
    )
    assert [c for c in columns_pairsam if c != "sam1" and c != "sam2"] == columns_pairs

    tmpdir.cleanup()
