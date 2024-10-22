# -*- coding: utf-8 -*-
import os
import sys

import pytest

import subprocess

testdir = os.path.dirname(os.path.realpath(__file__))


def test_mock_pysam_parse2_read():
    mock_sam_path = os.path.join(testdir, "data", "mock.parse2.sam")
    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "parse2",
                "-c",
                mock_chroms_path,
                "--add-pair-index",
                "--flip",
                "--report-position",
                "junction",
                "--report-orientation",
                "pair",
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
        print(l.split("SIMULATED:", 1)[1].split("\031", 1)[0].split("|"), id_counter)
        simulated_pair = (
            l.split("SIMULATED:", 1)[1]
            .split("\031", 1)[0]
            .split("|")[id_counter]
            .split(",")
        )
        print(assigned_pair)
        print(simulated_pair, prev_id)
        print()

        assert assigned_pair == simulated_pair


def test_mock_pysam_parse2_pair():
    mock_sam_path = os.path.join(testdir, "data", "mock.parse-all.sam")
    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "parse2",
                "-c",
                mock_chroms_path,
                "--add-pair-index",
                "--flip",
                "--report-position",
                "outer",
                "--report-orientation",
                "pair",
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
            l.split("SIMULATED:", 1)[1]
            .split("\031", 1)[0]
            .split("|")[id_counter]
            .split(",")
        )
        print(assigned_pair)
        print(simulated_pair, prev_id)
        print()

        assert assigned_pair == simulated_pair


def test_mock_pysam_parse2_single_end():
    
    """Testing single-end mode for parse2, no-flip mode.
    --report-position is outer (parse2 default)
    --report-orientation is pair (parse2 default)
    """

    mock_sam_path = os.path.join(testdir, "data", "mock.parse2-single-end.sam")
    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "parse2",
                "-c",
                mock_chroms_path,
                "--single-end",
                "--add-pair-index",
                "--no-flip",
                "--report-position",
                "outer",
                "--report-orientation",
                "pair",
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
        print(l.split("SIMULATED:", 1)[1].split("\031", 1)[0].split("|"), id_counter)
        simulated_pair = (
            l.split("SIMULATED:", 1)[1]
            .split("\031", 1)[0]
            .split("|")[id_counter]
            .split(",")
        )
        print(assigned_pair)
        print(simulated_pair, prev_id)
        print()

        assert assigned_pair == simulated_pair


def test_mock_pysam_parse2_single_end_expand():
    
    """Testing single-end mode for parse2, no-flip mode, with --expand.
    --report-position is outer (parse2 default)
    --report-orientation is pair (parse2 default)
    """

    mock_sam_path = os.path.join(testdir, "data", "mock.parse2-single-end.expand.sam")
    mock_chroms_path = os.path.join(testdir, "data", "mock.chrom.sizes")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "parse2",
                "-c",
                mock_chroms_path,
                "--single-end",
                "--expand",
                "--add-pair-index",
                "--no-flip",
                "--report-position",
                "outer",
                "--report-orientation",
                "pair",
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
        print(l.split("SIMULATED:", 1)[1].split("\031", 1)[0].split("|"), id_counter)
        simulated_pair = (
            l.split("SIMULATED:", 1)[1]
            .split("\031", 1)[0]
            .split("|")[id_counter]
            .split(",")
        )
        print(assigned_pair)
        print(simulated_pair, prev_id)
        print()

        assert assigned_pair == simulated_pair