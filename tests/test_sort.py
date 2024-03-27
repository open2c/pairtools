# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import pytest

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
        l.strip() for l in open(mock_pairsam_path, "r") if l.startswith("#")
    ]
    output_header = [l.strip() for l in result.split("\n") if l.startswith("#")]

    print(output_header)
    print(pairsam_header)
    for l in output_header:
        if not any([l in l2 for l2 in pairsam_header]):
            assert (
                l.startswith("#samheader: @PG")
                or l.startswith("#sorted")
                or l.startswith("#chromosomes")
            )

    pairsam_body = [
        l.strip()
        for l in open(mock_pairsam_path, "r")
        if not l.startswith("#") and l.strip()
    ]
    output_body = [
        l.strip() for l in result.split("\n") if not l.startswith("#") and l.strip()
    ]

    # check that all pairsam entries survived sorting:
    assert len(pairsam_body) == len(output_body)

    # check the sorting order of the output:
    prev_pair = None
    for l in output_body:
        cur_pair = l.split("\t")[1:8]
        if prev_pair is not None:
            assert cur_pair[0] >= prev_pair[0]
            if cur_pair[0] == prev_pair[0]:
                assert cur_pair[2] >= prev_pair[2]
                if cur_pair[2] == prev_pair[2]:
                    assert int(cur_pair[1]) >= int(prev_pair[1])
                    if int(cur_pair[1]) == int(prev_pair[1]):
                        assert int(cur_pair[3]) >= int(prev_pair[3])

        prev_pair = cur_pair
