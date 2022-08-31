# -*- coding: utf-8 -*-
import os
import sys

import pytest

import subprocess

testdir = os.path.dirname(os.path.realpath(__file__))


def test_restrict():
    """Restrict pairs file"""
    mock_pairs_path = os.path.join(testdir, "data", "mock.test-restr.pairs")
    mock_rfrag_path = os.path.join(testdir, "data", "mock.rsites.bed")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "restrict",
                "-f",
                mock_rfrag_path,
                mock_pairs_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    # check if the header got transferred correctly
    true_header = [l.strip() for l in open(mock_pairs_path, "r") if l.startswith("@")]
    output_header = [l.strip() for l in result.split("\n") if l.startswith("#")]
    for l in true_header:
        assert any([l in l2 for l2 in output_header])

    # check that the pairs got assigned properly
    cols = [x for x in output_header if x.startswith("#columns")][0].split(" ")[1:]

    COL_RFRAG1_TRUE = cols.index("rfrag_test1")
    COL_RFRAG2_TRUE = cols.index("rfrag_test2")
    COL_RFRAG1_OUTPUT = cols.index("rfrag1")
    COL_RFRAG2_OUTPUT = cols.index("rfrag2")

    for l in result.split("\n"):
        if l.startswith("#") or not l:
            continue

        line = l.split()
        assert line[COL_RFRAG1_TRUE] == line[COL_RFRAG1_OUTPUT]
        assert line[COL_RFRAG2_TRUE] == line[COL_RFRAG2_OUTPUT]
