# -*- coding: utf-8 -*-
import os
import subprocess
import sys
# import pytest

testdir = os.path.dirname(os.path.realpath(__file__))


def test_mock_pairsam():
    mock_pairsam_path = os.path.join(testdir, "data", "mock.pairsam")
    try:
        result = subprocess.check_output(
            ["python", "-m", "pairtools", "sort", mock_pairsam_path],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    # Check that the only changes strings are a @PG record of a SAM header,
    # the "#sorted" entry and chromosomes
    pairsam_header = [
        line.strip() for line in open(mock_pairsam_path, "r") if line.startswith("#")
    ]
    output_header = [line.strip() for line in result.split("\n") if line.startswith("#")]

    print(output_header)
    print(pairsam_header)
    for line in output_header:
        if not any([line in l2 for l2 in pairsam_header]):
            assert (
                line.startswith("#samheader: @PG")
                or line.startswith("#sorted")
                or line.startswith("#chromosomes")
            )

    pairsam_body = [
        line.strip()
        for line in open(mock_pairsam_path, "r")
        if not line.startswith("#") and line.strip()
    ]
    output_body = [
        line.strip() for line in result.split("\n") if not line.startswith("#") and line.strip()
    ]

    # check that all pairsam entries survived sorting:
    assert len(pairsam_body) == len(output_body)

    # check the sorting order of the output:
    prev_pair = None
    for line in output_body:
        cur_pair = line.split("\t")[1:8]
        if prev_pair is not None:
            assert cur_pair[0] >= prev_pair[0]
            if cur_pair[0] == prev_pair[0]:
                assert cur_pair[2] >= prev_pair[2]
                if cur_pair[2] == prev_pair[2]:
                    assert int(cur_pair[1]) >= int(prev_pair[1])
                    if int(cur_pair[1]) == int(prev_pair[1]):
                        assert int(cur_pair[3]) >= int(prev_pair[3])

        prev_pair = cur_pair
