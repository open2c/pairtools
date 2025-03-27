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

    # Check that all pairsam entries survived merging
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

    # Check the sorting order of the output
    keys = [
        (pair[1], int(pair[2]), pair[3], int(pair[4]), pair[7])
        for pair in (l.split("\t") for l in output_body)
    ]
    assert all(keys[i] <= keys[i + 1] for i in range(len(keys) - 1)), "Output not sorted correctly"

    # Check that the header is preserved
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

    # Check the headers
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

    assert len(pairsam_header_1) + 1 == len(output_header)  # +1 for the PG line

    tmpdir.cleanup()