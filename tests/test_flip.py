# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import pytest

testdir = os.path.dirname(os.path.realpath(__file__))
mock_pairs_path = os.path.join(testdir, "data", "mock.4flip.pairs")
mock_chromsizes_path = os.path.join(testdir, "data", "mock.chrom.sizes")


def test_flip():
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "flip",
                mock_pairs_path,
                "-c",
                mock_chromsizes_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    orig_pairs = [
        l.strip().split("\t")
        for l in open(mock_pairs_path, "r")
        if not l.startswith("#") and l.strip()
    ]
    flipped_pairs = [
        l.strip().split("\t")
        for l in result.split("\n")
        if not l.startswith("#") and l.strip()
    ]

    chrom_enum = {"!": 0, "chr1": 1, "chr2": 2, "chrU": 3, "chrU1": 4}
    # chrU stands for unannotated chromosome, which has less priority than annotated ones
    # chrU1 is another unannotated chromosome, which should go lexigographically after chrU
    for orig_pair, flipped_pair in zip(orig_pairs, flipped_pairs):
        has_correct_order = (chrom_enum[orig_pair[1]], int(orig_pair[2])) <= (
            chrom_enum[orig_pair[3]],
            int(orig_pair[4]),
        )
        if has_correct_order:
            assert all([c1 == c2 for c1, c2 in zip(orig_pair, flipped_pair)])
        if not has_correct_order:
            assert orig_pair[1] == flipped_pair[3]
            assert orig_pair[2] == flipped_pair[4]
            assert orig_pair[3] == flipped_pair[1]
            assert orig_pair[4] == flipped_pair[2]
            assert orig_pair[5] == flipped_pair[6]
            assert orig_pair[6] == flipped_pair[5]
            assert orig_pair[7] == flipped_pair[7][::-1]
