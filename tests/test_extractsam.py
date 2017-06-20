# -*- coding: utf-8 -*-
import os
import sys
import subprocess
from nose.tools import assert_raises, with_setup
import tempfile

testdir = os.path.dirname(os.path.realpath(__file__))
mock_pairsam_path = os.path.join(testdir, 'data', 'mock.pairsam')

tmpdir = tempfile.TemporaryDirectory()
tmpdir_name = tmpdir.name
pairs_path = os.path.join(tmpdir_name, 'out.pairs')
sam_path   = os.path.join(tmpdir_name, 'out.sam')

def setup_func():
    try:
        subprocess.check_output(
            ['python',
             '-m',
             'pairtools',
             'extractsam',
             mock_pairsam_path,
             '--output-pairs',
             pairs_path,
             '--output-sam',
             sam_path,
             ],
            )
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

def teardown_func():
    tmpdir.cleanup()

@with_setup(setup_func, teardown_func)
def test_extractsam():

    pairsam_lines = [l.strip() for l in open(mock_pairsam_path, 'r') 
                     if l.strip()]
    pairs_lines = [l.strip() for l in open(pairs_path, 'r')
                   if l.strip()]
    sam_lines = [l.strip() for l in open(sam_path, 'r')
                if l.strip()]

    # check that all entries survived splitting:
    n_pairsam = len([l for l in pairsam_lines if not l.startswith('#')])
    n_pairs = len([l for l in pairs_lines if not l.startswith('#')])
    n_sam =  len([l for l in sam_lines if not l.startswith('@')])  // 2

    assert n_pairsam == n_pairs
    assert n_pairsam == n_sam

    # check that the header survived splitting:
    pairsam_header = [l.strip() for l in open(mock_pairsam_path, 'r') 
                     if l.strip() and l.startswith('#')]
    pairs_header = [l.strip() for l in open(pairs_path, 'r')
                    if l.strip() and l.startswith('#')]
    sam_header = [l.strip() for l in open(sam_path, 'r')
                  if l.strip() and l.startswith('@')]
    assert all(
            any(l in l2 for l2 in pairsam_header)
            for l in sam_header if not l.startswith('@PG'))


