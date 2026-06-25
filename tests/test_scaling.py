# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import pytest
import pandas as pd
import io

testdir = os.path.dirname(os.path.realpath(__file__))


def test_scaling():
    mock_pairsam_path = os.path.join(testdir, "data", "mock.pairsam")
    try:
        result = subprocess.check_output(
            ["python", "-m", "pairtools", "scaling", mock_pairsam_path],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    output = pd.read_csv(io.StringIO(result), sep="\t", header=0)

    assert (
        output["n_pairs"].sum() == 9
    )  # double unmapped pairs are currently ignored by lib.scaling


def test_scaling_uses_header_chromsizes():
    # Without an explicit view, the implicit per-chromosome view must use the chromsizes
    # from the .pairs header (chr* = 100 in mock.pairsam); otherwise region ends stay at
    # -1 and n_bp2 is left empty.
    from pairtools.lib.scaling import compute_scaling

    mock_pairsam_path = os.path.join(testdir, "data", "mock.pairsam")
    cis, _trans = compute_scaling(mock_pairsam_path, regions=None)

    assert (cis["end1"] == 100).all()
    assert (cis["n_bp2"] > 0).any()
