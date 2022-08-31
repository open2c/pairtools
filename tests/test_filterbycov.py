# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import pytest
import tempfile

testdir = os.path.dirname(os.path.realpath(__file__))

mock_pairs_path_filterbycov = os.path.join(testdir, "data", "mock.4filterbycov.pairs")

tmpdir = tempfile.TemporaryDirectory()
tmpdir_name = tmpdir.name


params = [
    {"max_dist": 0, "max_cov": 3},
    {"max_dist": 0, "max_cov": 2},
    {"max_dist": 1, "max_cov": 1},
]

for p in params:
    p["lowcov_path"] = os.path.join(
        tmpdir_name, "lowcov.{}.{}.pairs".format(p["max_dist"], p["max_cov"])
    )
    p["highcov_path"] = os.path.join(
        tmpdir_name, "highcov.{}.{}.pairs".format(p["max_dist"], p["max_cov"])
    )
    p["unmapped_path"] = os.path.join(
        tmpdir_name, "unmapped.{}.{}.pairs".format(p["max_dist"], p["max_cov"])
    )


@pytest.fixture
def setup_filterbycov():
    try:
        for p in params:
            subprocess.check_output(
                [
                    "python",
                    "-m",
                    "pairtools",
                    "filterbycov",
                    mock_pairs_path_filterbycov,
                    "--output",
                    p["lowcov_path"],
                    "--output-highcov",
                    p["highcov_path"],
                    "--output-unmapped",
                    p["unmapped_path"],
                    "--max-dist",
                    str(p["max_dist"]),
                    "--max-cov",
                    str(p["max_cov"]),
                ]
            )
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e


def test_mock_pairs(setup_filterbycov):

    all_pairs = [
        l.strip().split("\t")
        for l in open(mock_pairs_path_filterbycov, "r")
        if not l.startswith("#") and l.strip()
    ]
    for p in params:

        lowcov_pairs = [
            l.strip().split("\t")
            for l in open(p["lowcov_path"], "r")
            if not l.startswith("#") and l.strip()
        ]
        highcov_pairs = [
            l.strip().split("\t")
            for l in open(p["highcov_path"], "r")
            if not l.startswith("#") and l.strip()
        ]
        unmapped_pairs = [
            l.strip().split("\t")
            for l in open(p["unmapped_path"], "r")
            if not l.startswith("#") and l.strip()
        ]

        # check that at least a few pairs remained in deduped and dup files
        # assert len(lowcov_pairs) > 0
        # assert len(highcov_pairs) > 0
        # assert len(unmapped_pairs) > 0

        # check that all pairs entries survived deduping:

        assert len(lowcov_pairs) + len(unmapped_pairs) + len(highcov_pairs) == len(
            all_pairs
        )

        assert all([(pair[1] != "!" and pair[3] != "!") for pair in lowcov_pairs])
        assert all([(pair[1] != "!" and pair[3] != "!") for pair in highcov_pairs])
        assert all([(pair[1] == "!" or pair[3] == "!") for pair in unmapped_pairs])

        def update_coverage(coverage, chrom, pos, max_dist):
            if chrom == "!":
                return

            coverage[chrom] = coverage.get(chrom, {})
            for i in range(max(0, pos - max_dist), pos + max_dist + 1):
                coverage[chrom][i] = coverage[chrom].get(i, 0) + 1

        coverage = {}
        for pair in all_pairs:
            update_coverage(coverage, pair[1], int(pair[2]), p["max_dist"])
            update_coverage(coverage, pair[3], int(pair[4]), p["max_dist"])

        for pair in lowcov_pairs:
            # print (p['max_cov'],p['max_dist'])
            # print (pair, coverage[pair[1]][int(pair[2])])
            # print (pair, coverage[pair[3]][int(pair[4])])
            assert coverage[pair[1]][int(pair[2])] <= p["max_cov"]
            assert coverage[pair[3]][int(pair[4])] <= p["max_cov"]

        for pair in highcov_pairs:
            # print (p['max_cov'],p['max_dist'])
            # print (pair, coverage[pair[1]][int(pair[2])])
            # print (pair, coverage[pair[3]][int(pair[4])])
            assert (coverage[pair[1]][int(pair[2])] > p["max_cov"]) or (
                coverage[pair[3]][int(pair[4])] > p["max_cov"]
            )

    tmpdir.cleanup()
