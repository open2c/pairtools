# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import pytest


import tempfile

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
        l.strip()
        for l in open(mock_pairsam_path_1, "r")
        if not l.startswith("#") and l.strip()
    ]
    pairsam_body_2 = [
        l.strip()
        for l in open(mock_pairsam_path_2, "r")
        if not l.startswith("#") and l.strip()
    ]
    output_body = [
        l.strip() for l in result.split("\n") if not l.startswith("#") and l.strip()
    ]
    assert len(pairsam_body_1) + len(pairsam_body_2) == len(output_body)

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
        l.strip()
        for l in open(mock_sorted_pairsam_path_1, "r")
        if l.startswith("#") and l.strip()
    ]
    pairsam_header_2 = [
        l.strip()
        for l in open(mock_sorted_pairsam_path_2, "r")
        if l.startswith("#") and l.strip()
    ]
    output_header = [
        l.strip() for l in result.split("\n") if l.startswith("#") and l.strip()
    ]

    assert len(pairsam_header_1) + 1 == len(output_header)

    tmpdir.cleanup()
