# -*- coding: utf-8 -*-
import os
import subprocess
import sys

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
        print(line.split("SIMULATED:", 1)[1].split("\031", 1)[0].split("|"), id_counter)
        simulated_pair = (
            line.split("SIMULATED:", 1)[1]
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
            line.split("SIMULATED:", 1)[1]
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
        print(line.split("SIMULATED:", 1)[1].split("\031", 1)[0].split("|"), id_counter)
        simulated_pair = (
            line.split("SIMULATED:", 1)[1]
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
        print(line.split("SIMULATED:", 1)[1].split("\031", 1)[0].split("|"), id_counter)
        simulated_pair = (
            line.split("SIMULATED:", 1)[1]
            .split("\031", 1)[0]
            .split("|")[id_counter]
            .split(",")
        )
        print(assigned_pair)
        print(simulated_pair, prev_id)
        print()

        assert assigned_pair == simulated_pair
