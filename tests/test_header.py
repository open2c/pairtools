# -*- coding: utf-8 -*-
import os
import sys

from nose.tools import assert_raises

import subprocess
from pairtools import _pairsam_format

testdir = os.path.dirname(os.path.realpath(__file__))


def test_generate():
    """Test generation of the header.
    Example run:
    pairtools header generate tests/data/mock.pairsam \
    --chroms-path tests/data/mock.chrom.sizes --pairsam \
    --sam-path tests/data/mock.sam
    """

    mock_sam_path = os.path.join(testdir, "data", "mock.sam")
    mock_pairs_path = os.path.join(testdir, "data", "mock.pairsam")
    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "header",
                "generate",
                "--chroms-path",
                mock_chroms_path,
                "--sam-path",
                mock_sam_path,
                "--pairsam",
                mock_pairs_path,
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


def test_remove():
    """Test removal of columns from the file
    Example run:
    pairtools header remove-columns tests/data/mock.pairsam -c sam1,sam2
    """

    mock_pairs_path = os.path.join(testdir, "data", "mock.pairsam")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "header",
                "remove-columns",
                "-c",
                "sam1,sam2",
                mock_pairs_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    # check if the columns are removed properly:
    pairsam_header = [l.strip() for l in result.split("\n") if l.startswith("#")]
    for l in pairsam_header:
        if l.startswith("#columns:"):
            line = l.strip()
            assert (
                line
                == "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type"
            )

    # check that the pairs got assigned properly
    for l in result.split("\n"):
        if l.startswith("#") or not l:
            continue

        assert len(l.split(_pairsam_format.PAIRSAM_SEP)) == 8
