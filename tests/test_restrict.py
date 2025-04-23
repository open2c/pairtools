# -*- coding: utf-8 -*-
import os
import subprocess
import sys

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
    true_header = [line.strip() for line in open(mock_pairs_path, "r") if line.startswith("@")]
    output_header = [line.strip() for line in result.split("\n") if line.startswith("#")]
    for header_line in true_header:
        assert any([header_line in output_line for output_line in output_header])

    # check that the pairs got assigned properly
    cols = [x for x in output_header if x.startswith("#columns")][0].split(" ")[1:]

    COL_RFRAG1_TRUE = cols.index("rfrag_test1")
    COL_RFRAG2_TRUE = cols.index("rfrag_test2")
    COL_RFRAG1_OUTPUT = cols.index("rfrag1")
    COL_RFRAG2_OUTPUT = cols.index("rfrag2")

    for line in result.split("\n"):
        if line.startswith("#") or not line:
            continue

        line_data = line.split()
        assert line_data[COL_RFRAG1_TRUE] == line_data[COL_RFRAG1_OUTPUT]
        assert line_data[COL_RFRAG2_TRUE] == line_data[COL_RFRAG2_OUTPUT]
