# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import click
from nose.tools import assert_raises, with_setup
import tempfile

testdir = os.path.dirname(os.path.realpath(__file__))
mock_pairsam_path = os.path.join(testdir, 'data', 'mock.4dedup.pairsam')

tmpdir = tempfile.TemporaryDirectory()
tmpdir_name = tmpdir.name
dedup_path = os.path.join(tmpdir_name, 'dedup.pairsam')
dups_path = os.path.join(tmpdir_name, 'dups.pairsam')

max_mismatch = 3
def setup_func():
    try:
        subprocess.check_output(
            ['python',
             '-m',
             'pairsamtools',
             'dedup',
             '--input',
             mock_pairsam_path,
             '--output',
             dedup_path,
             '--output-dups',
             dups_path,
             '--max-mismatch',
             str(max_mismatch)
             ],
            )
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

def teardown_func():
    tmpdir.cleanup()

@with_setup(setup_func, teardown_func)
def test_mock_pairsam():

    pairsam_body = [l.strip() for l in open(mock_pairsam_path, 'r') 
                      if not l.startswith('#') and l.strip()]
    output_body  = [l.strip() for l in open(dedup_path, 'r')
                    if not l.startswith('#') and l.strip()]
    output_dups_body  = [l.strip() for l in open(dups_path, 'r')
                         if not l.startswith('#') and l.strip()]

    # check that at least a few pairs remained in deduped and dup files
    assert len(output_body) > 0
    assert len(output_dups_body) > 0

    # check that all pairsam entries survived deduping:
    assert len(output_body) + len(output_dups_body) == len(pairsam_body)


    dedup_pairs = [l.split('\v') for l in output_body]
    dup_pairs   = [l.split('\v') for l in output_dups_body]

    def pairs_overlap(pair1, pair2, max_mismatch):
        overlap = (
            (pair1[1] == pair2[1])
            and (pair1[2] == pair2[2])
            and (pair1[5] == pair2[5])
            and (pair1[6] == pair2[6])
            and (abs(int(pair1[3]) - int(pair2[3])) <= max_mismatch)
            and (abs(int(pair1[4]) - int(pair2[4])) <= max_mismatch)
            )
        return overlap

    # check that deduped pairs do not overlap
    assert all([not pairs_overlap(pair1, pair2, 3)
                for i, pair1 in enumerate(dedup_pairs)
                for j, pair2 in enumerate(dedup_pairs)
                if i != j])

    # check that the removed duplicates overlap with at least one of the 
    # deduplicated entries
    assert all([
            any([pairs_overlap(pair1, pair2, 3)
                 for pair2 in dedup_pairs])
                for pair1 in dup_pairs
            ])

