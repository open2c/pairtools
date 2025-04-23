# -*- coding: utf-8 -*-
import os
import subprocess
import sys

# import pytest

testdir = os.path.dirname(os.path.realpath(__file__))


def test_python_version():
    assert sys.version_info[0] == 3, "Use Python 3!"


def test_mock_pysam():
    """Parse non-chimeric alignments with walks-policy mask with pysam backend."""
    mock_sam_path = os.path.join(testdir, "data", "mock.sam")
    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "parse",
                "--walks-policy",
                "mask",
                "-c",
                mock_chroms_path,
                mock_sam_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    # check if the header got transferred correctly
    sam_header = [line.strip() for line in open(mock_sam_path, "r") if line.startswith("@")]
    pairsam_header = [line.strip() for line in result.split("\n") if line.startswith("#")]
    for header_line in sam_header:
        assert any([header_line in pairsam_line for pairsam_line in pairsam_header])

    # check that the pairs got assigned properly
    for line in result.split("\n"):
        if line.startswith("#") or not line:
            continue

        print(line)
        assigned_pair = line.split("\t")[1:8]
        simulated_pair = line.split("CT:Z:SIMULATED:", 1)[1].split("\031", 1)[0].split(",")
        print(assigned_pair)
        print(simulated_pair)
        print()

        assert assigned_pair == simulated_pair


def test_mock_pysam_parse_all():
    """Parse all alignment in each read with walks-policy all and pysam backend."""
    mock_sam_path = os.path.join(testdir, "data", "mock.parse-all.sam")
    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "parse",
                "--walks-policy",
                "all",
                "-c",
                mock_chroms_path,
                "--add-pair-index",
                mock_sam_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    # check if the header got transferred correctly
    sam_header = [line.strip() for line in open(mock_sam_path, "r") if line.startswith("@")]
    pairsam_header = [line.strip() for line in result.split("\n") if line.startswith("#")]
    for header_line in sam_header:
        assert any([header_line in pairsam_line for pairsam_line in pairsam_header])

    # check that the pairs got assigned properly
    id_counter = 0
    prev_id = ""
    for line in result.split("\n"):
        if line.startswith("#") or not line:
            continue

        if prev_id == line.split("\t")[0]:
            id_counter += 1
        else:
            id_counter = 0
        prev_id = line.split("\t")[0]

        assigned_pair = line.split("\t")[1:8] + line.split("\t")[-2:]
        simulated_pair = (
            line.split("CT:Z:SIMULATED:", 1)[1]
            .split("\031", 1)[0]
            .split("|")[id_counter]
            .split(",")
        )
        print(assigned_pair)
        print(simulated_pair, prev_id)
        print()

        assert assigned_pair == simulated_pair
