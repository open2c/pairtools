# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import pytest
import tempfile

testdir = os.path.dirname(os.path.realpath(__file__))


def test_mock_pairsam(setup_dedup):
    mock_pairsam_path, dedup_path, unmapped_path, dups_path, \
    dedup_max_path, unmapped_max_path, dups_max_path, \
    dedup_markdups_path, unmapped_markdups_path, dups_markdups_path, max_mismatch, tmpdir = setup_dedup

    pairsam_pairs = [
        l.strip().split("\t")
        for l in open(mock_pairsam_path, "r")
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