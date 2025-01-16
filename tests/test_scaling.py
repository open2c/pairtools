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

    assert output["n_pairs"].sum() == 7  # unmapped pairs are ignored by lib.scaling
