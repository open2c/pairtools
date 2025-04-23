# -*- coding: utf-8 -*-
import os
import subprocess
import sys
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


def test_custom_column_warning(tmpdir):
    """Test that a warning is emitted when sorting with a custom column not in pairsam format."""
    # Create a temporary .pairsam file with custom_col in the header
    mock_pairsam_path = os.path.join(tmpdir, "test.pairsam")
    output_path = os.path.join(tmpdir, "sorted_output.pairsam")

    # Write a minimal .pairsam file
    with open(mock_pairsam_path, "w") as f:
        f.write("## pairs format v1.0\n")
        f.write("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 pair_type custom_col\n")
        f.write("read1\tchr1\t100\tchr2\t200\t+\t-\tUU\t42\n")
        f.write("read2\tchr2\t150\tchr1\t250\t-\t+\tUU\t99\n")

    # Run sort command with a custom column
    cmd = [
        "python",
        "-m",
        "pairtools",
        "sort",
        mock_pairsam_path,
        "--output",
        output_path,
        "--extra-col",
        "custom_col",
    ]

    # Capture stderr to check for warning
    process = subprocess.Popen(
        cmd,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    stdout, stderr = process.communicate()

    # Check that the command completed successfully
    assert process.returncode == 0, f"Command failed: {stderr}"

    # Verify warning was emitted
    assert "Column 'custom_col' not found in pairsam format definitions" in stderr
    assert "Assuming string type for sorting" in stderr

    # Verify output file exists and has content
    assert os.path.exists(output_path)
    with open(output_path, "r") as f:
        output_lines = f.readlines()
    assert len(output_lines) > 0

    # Check that the output is sorted (basic check on header and body)
    output_header = [line.strip() for line in output_lines if line.startswith("#")]
    output_body = [
        line.strip() for line in output_lines if not line.startswith("#") and line.strip()
    ]
    assert any(line.startswith("#sorted") for line in output_header)
    assert len(output_body) > 0
