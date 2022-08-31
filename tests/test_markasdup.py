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
            ["python", "-m", "pairtools", "markasdup", mock_pairsam_path],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

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

    # check that all pairtypes got changed to DD
    for l in output_body:
        assert l.split("\t")[7] == "DD"
