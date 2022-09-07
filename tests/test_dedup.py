# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import pytest
import tempfile

testdir = os.path.dirname(os.path.realpath(__file__))

tmpdir = tempfile.TemporaryDirectory()
tmpdir_name = tmpdir.name

mock_pairsam_path_dedup = os.path.join(testdir, "data", "mock.4dedup.pairsam")

dedup_path = os.path.join(tmpdir_name, "dedup.pairsam")
unmapped_path = os.path.join(tmpdir_name, "unmapped.pairsam")
dups_path = os.path.join(tmpdir_name, "dups.pairsam")

dedup_max_path = os.path.join(tmpdir_name, "dedup_max.pairsam")
unmapped_max_path = os.path.join(tmpdir_name, "unmapped_max.pairsam")
dups_max_path = os.path.join(tmpdir_name, "dups_max.pairsam")

dedup_markdups_path = os.path.join(tmpdir_name, "dedup.markdups.pairsam")
unmapped_markdups_path = os.path.join(tmpdir_name, "unmapped.markdups.pairsam")
dups_markdups_path = os.path.join(tmpdir_name, "dups.markdups.pairsam")

max_mismatch = 1


@pytest.fixture
def setup_dedup():
    try:
        subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "dedup",
                mock_pairsam_path_dedup,
                "--output",
                dedup_path,
                "--output-dups",
                dups_path,
                "--output-unmapped",
                unmapped_path,
                "--max-mismatch",
                str(max_mismatch),
            ],
        )
        subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "dedup",
                mock_pairsam_path_dedup,
                "--output",
                dedup_max_path,
                "--output-dups",
                dups_max_path,
                "--output-unmapped",
                unmapped_max_path,
                "--max-mismatch",
                str(max_mismatch),
                "--method",
                "max",
            ],
        )
        subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "dedup",
                mock_pairsam_path_dedup,
                "--mark-dups",
                "--output",
                dedup_markdups_path,
                "--output-dups",
                dups_markdups_path,
                "--output-unmapped",
                unmapped_markdups_path,
                "--max-mismatch",
                str(max_mismatch),
            ],
        )
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e


def test_mock_pairsam(setup_dedup):

    pairsam_pairs = [
        l.strip().split("\t")
        for l in open(mock_pairsam_path_dedup, "r")
        if not l.startswith("#") and l.strip()
    ]
    for (ddp, up, dp) in [
        (dedup_path, unmapped_path, dups_path),
        (dedup_max_path, unmapped_max_path, dups_max_path),
        (dedup_markdups_path, unmapped_markdups_path, dups_markdups_path),
    ]:

        dedup_pairs = [
            l.strip().split("\t")
            for l in open(ddp, "r")
            if not l.startswith("#") and l.strip()
        ]
        unmapped_pairs = [
            l.strip().split("\t")
            for l in open(up, "r")
            if not l.startswith("#") and l.strip()
        ]
        dup_pairs = [
            l.strip().split("\t")
            for l in open(dp, "r")
            if not l.startswith("#") and l.strip()
        ]

        # check that at least a few pairs remained in deduped and dup files
        assert len(dedup_pairs) > 0
        assert len(dup_pairs) > 0
        assert len(unmapped_pairs) > 0
        import pandas as pd

        # check that all pairsam entries survived deduping:

        assert len(dedup_pairs) + len(unmapped_pairs) + len(dup_pairs) == len(
            pairsam_pairs
        )

        def pairs_overlap(pair1, pair2, max_mismatch):
            overlap = (
                (pair1[1] == pair2[1])
                and (pair1[3] == pair2[3])
                and (pair1[5] == pair2[5])
                and (pair1[6] == pair2[6])
                and (abs(int(pair1[2]) - int(pair2[2])) <= max_mismatch)
                and (abs(int(pair1[4]) - int(pair2[4])) <= max_mismatch)
            )
            return overlap

        # check that deduped pairs do not overlap
        assert all(
            [
                not pairs_overlap(pair1, pair2, max_mismatch)
                for i, pair1 in enumerate(dedup_pairs)
                for j, pair2 in enumerate(dedup_pairs)
                if i != j
            ]
        )

        # check that the removed duplicates overlap with at least one of the
        # deduplicated entries
        assert all(
            [
                any([pairs_overlap(pair1, pair2, 3) for pair2 in dedup_pairs])
                for pair1 in dup_pairs
            ]
        )

    tmpdir.cleanup()
