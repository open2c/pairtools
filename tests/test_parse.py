# -*- coding: utf-8 -*-
import os
import sys

import pytest

import subprocess

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
    sam_header = [l.strip() for l in open(mock_sam_path, "r") if l.startswith("@")]
    pairsam_header = [l.strip() for l in result.split("\n") if l.startswith("#")]
    for l in sam_header:
        assert any([l in l2 for l2 in pairsam_header])

    # check that the pairs got assigned properly
    for l in result.split("\n"):
        if l.startswith("#") or not l:
            continue

        print(l)
        assigned_pair = l.split("\t")[1:8]
        simulated_pair = l.split("CT:Z:SIMULATED:", 1)[1].split("\031", 1)[0].split(",")
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
    sam_header = [l.strip() for l in open(mock_sam_path, "r") if l.startswith("@")]
    pairsam_header = [l.strip() for l in result.split("\n") if l.startswith("#")]
    for l in sam_header:
        assert any([l in l2 for l2 in pairsam_header])

    # check that the pairs got assigned properly
    id_counter = 0
    prev_id = ""
    for l in result.split("\n"):
        if l.startswith("#") or not l:
            continue

        if prev_id == l.split("\t")[0]:
            id_counter += 1
        else:
            id_counter = 0
        prev_id = l.split("\t")[0]

        assigned_pair = l.split("\t")[1:8] + l.split("\t")[-2:]
        simulated_pair = (
            l.split("CT:Z:SIMULATED:", 1)[1]
            .split("\031", 1)[0]
            .split("|")[id_counter]
            .split(",")
        )
        print(assigned_pair)
        print(simulated_pair, prev_id)
        print()

        assert assigned_pair == simulated_pair
