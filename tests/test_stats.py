# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import numpy as np
import yaml

import pytest


testdir = os.path.dirname(os.path.realpath(__file__))


def test_mock_pairsam():
    mock_pairsam_path = os.path.join(testdir, "data", "mock.4stats.pairs")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "stats",
                "--yaml",
                "--single-mapped-by-side",
                mock_pairsam_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    stats = yaml.safe_load(result)

    # for k in stats["no_filter"]:
    #     try:
    #         stats["no_filter"][k] = int(stats["no_filter"][k])
    #     except (ValueError, TypeError):
    #         stats["no_filter"][k] = float(stats["no_filter"][k])

    assert stats["no_filter"]["total"] == 9
    assert stats["no_filter"]["total_single_sided_mapped"] == 2
    assert stats["no_filter"]["total_left_only_mapped"] == 0
    assert stats["no_filter"]["total_right_only_mapped"] == 2
    assert stats["no_filter"]["total_mapped"] == 6
    assert stats["no_filter"]["total_dups"] == 1
    assert stats["no_filter"]["cis"] == 3
    assert stats["no_filter"]["trans"] == 2
    assert stats["no_filter"]["pair_types"]["UU"] == 4
    assert stats["no_filter"]["pair_types"]["NU"] == 1
    assert stats["no_filter"]["pair_types"]["WW"] == 1
    assert stats["no_filter"]["pair_types"]["UR"] == 1
    assert stats["no_filter"]["pair_types"]["MU"] == 1
    assert stats["no_filter"]["pair_types"]["DD"] == 1
    assert stats["no_filter"]["chrom_freq"]["chr1/chr2"] == 1
    assert stats["no_filter"]["chrom_freq"]["chr1/chr1"] == 3
    assert stats["no_filter"]["chrom_freq"]["chr2/chr3"] == 1
    for orientation in ("++", "+-", "-+", "--"):
        s = stats["no_filter"]["dist_freq"][orientation]
        for k, val in s.items():
            if orientation == "++" and k in [1, 2, 42]:
                assert s[k] == 1
            else:
                assert s[k] == 0

    assert stats["no_filter"]["summary"]["frac_cis"] == 0.6
    assert stats["no_filter"]["summary"]["frac_cis_1kb+"] == 0
    assert stats["no_filter"]["summary"]["frac_cis_2kb+"] == 0
    assert stats["no_filter"]["summary"]["frac_cis_4kb+"] == 0
    assert stats["no_filter"]["summary"]["frac_cis_10kb+"] == 0
    assert stats["no_filter"]["summary"]["frac_cis_20kb+"] == 0
    assert stats["no_filter"]["summary"]["frac_cis_40kb+"] == 0
    assert np.isclose(stats["no_filter"]["summary"]["frac_dups"], 1 / 6)


def test_merge_stats():
    mock_pairsam_path = os.path.join(testdir, "data", "mock.4stats.pairs")
    try:
        subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "stats",
                "--with-chromsizes",
                mock_pairsam_path,
                "--output",
                "mock.stats",
            ],
        )

        subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "stats",
                "--no-chromsizes",
                mock_pairsam_path,
                "--output",
                "mock.no_chromsizes.stats",
            ],
        )
        subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "stats",
                "mock.stats",
                "mock.stats",
                "--merge",
                "--output",
                "mock.merged_chromsizes.stats",
            ],
        )
        subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "stats",
                "mock.stats",
                "mock.no_chromsizes.stats",
                "--merge",
                "--output",
                "mock.merged_mixed.stats",
            ],
        )
        subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "stats",
                "mock.no_chromsizes.stats",
                "mock.no_chromsizes.stats",
                "--merge",
                "--output",
                "mock.merged_no_chromsizes.stats",
            ],
        )
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e


from pairtools.lib.stats import PairCounter


@pytest.fixture
def pair_counter():
    counter = PairCounter(filters={"f1": "filter1", "f2": "filter2"})
    counter._dist_bins = np.array([1, 1000, 10000, 100000, 1000000])
    # Populate the counter with some sample data
    counter._stat["f1"]["dist_freq"] = {
        "++": {1: 80, 1000: 80, 10000: 91, 100000: 95},
        "--": {1: 100, 1000: 100, 10000: 100, 100000: 100},
        "-+": {1: 100, 1000: 100, 10000: 100, 100000: 100},
        "+-": {1: 120, 1000: 120, 10000: 109, 100000: 105},
    }

    counter._stat["f2"]["dist_freq"] = {
        "++": {1: 200, 1000: 180, 10000: 160, 100000: 140},
        "--": {1: 220, 1000: 190, 10000: 170, 100000: 150},
        "-+": {1: 210, 1000: 185, 10000: 165, 100000: 145},
        "+-": {1: 230, 1000: 195, 10000: 175, 100000: 155},
    }

    return counter


def test_find_dist_freq_convergence_distance(pair_counter):
    result = pair_counter.find_dist_freq_convergence_distance(0.1)

    assert "f1" in result
    assert "f2" in result

    f1_result = result["f1"]
    assert "convergence_dist" in f1_result
    assert "strands_w_max_convergence_dist" in f1_result
    assert "convergence_rel_diff_threshold" in f1_result
    assert "n_cis_pairs_below_convergence_dist" in f1_result
    assert "n_cis_pairs_below_convergence_dist_all_strands" in f1_result
    assert "n_cis_pairs_above_convergence_dist_all_strands" in f1_result
    assert "frac_cis_in_cis_below_convergence_dist" in f1_result
    assert "frac_cis_in_cis_below_convergence_dist_all_strands" in f1_result
    assert "frac_cis_in_cis_above_convergence_dist_all_strands" in f1_result
    assert "frac_total_mapped_in_cis_below_convergence_dist" in f1_result
    assert "frac_total_mapped_in_cis_below_convergence_dist_all_strands" in f1_result
    assert "frac_total_mapped_in_cis_above_convergence_dist_all_strands" in f1_result

    assert f1_result["convergence_rel_diff_threshold"] == 0.1
    assert f1_result["convergence_dist"] == 10000
    assert f1_result["strands_w_max_convergence_dist"] == "++"

    # f2_result = result["f2"]
    # assert "convergence_dist" in f2_result
    # assert "strands_w_max_convergence_dist" in f2_result
    # assert "convergence_rel_diff_threshold" in f2_result
    # Add more assertions for f2_result as needed
