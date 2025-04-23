# -*- coding: utf-8 -*-
import os
import subprocess
import sys
import tempfile

import pytest

testdir = os.path.dirname(os.path.realpath(__file__))

tmpdir = tempfile.TemporaryDirectory()
tmpdir_name = tmpdir.name
mock_pairsam_path_1 = os.path.join(testdir, "data", "mock.pairsam")
mock_pairsam_path_2 = os.path.join(testdir, "data", "mock.2.pairsam")
mock_sorted_pairsam_path_1 = os.path.join(tmpdir_name, "1.pairsam")
mock_sorted_pairsam_path_2 = os.path.join(tmpdir_name, "2.pairsam")


@pytest.fixture
def setup_sort_two():
    try:
        subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "sort",
                mock_pairsam_path_1,
                "--output",
                mock_sorted_pairsam_path_1,
            ],
        )

        subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "sort",
                mock_pairsam_path_2,
                "--output",
                mock_sorted_pairsam_path_2,
            ],
        )
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e


def test_mock_pairsam(setup_sort_two):

    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "merge",
                mock_sorted_pairsam_path_1,
                mock_sorted_pairsam_path_2,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    # check that all pairsam entries survived sorting:
    pairsam_body_1 = [
        line.strip()
        for line in open(mock_pairsam_path_1, "r")
        if not line.startswith("#") and line.strip()
    ]
    pairsam_body_2 = [
        line.strip()
        for line in open(mock_pairsam_path_2, "r")
        if not line.startswith("#") and line.strip()
    ]
    output_body = [
        line.strip() for line in result.split("\n") if not line.startswith("#") and line.strip()
    ]
    assert len(pairsam_body_1) + len(pairsam_body_2) == len(output_body)

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

    # Check that the header is preserved:
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "merge",
                "--keep-first-header",
                mock_sorted_pairsam_path_1,
                mock_sorted_pairsam_path_2,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    # check the headers:
    pairsam_header_1 = [
        line.strip()
        for line in open(mock_sorted_pairsam_path_1, "r")
        if line.startswith("#") and line.strip()
    ]
    output_header = [
        line.strip() for line in result.split("\n") if line.startswith("#") and line.strip()
    ]

    assert len(pairsam_header_1) + 1 == len(output_header)

    tmpdir.cleanup()
