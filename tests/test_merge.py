# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import click
from nose.tools import assert_raises, with_setup
import tempfile

testdir = os.path.dirname(os.path.realpath(__file__))

tmpdir = tempfile.TemporaryDirectory()
tmpdir_name = tmpdir.name
mock_pairsam_path_1 = os.path.join(testdir, 'data', 'mock.pairsam')
mock_pairsam_path_2 = os.path.join(testdir, 'data', 'mock.2.pairsam')
mock_sorted_pairsam_path_1 = os.path.join(tmpdir_name, '1.pairsam')
mock_sorted_pairsam_path_2 = os.path.join(tmpdir_name, '2.pairsam')

def setup_func():
    subprocess.check_output(
        ['python',
         '-m',
         'pairsamtools',
         'sort',
         '--input',
         mock_pairsam_path_1,
         '--output',
         mock_sorted_pairsam_path_1
         ],
        )

    subprocess.check_output(
        ['python',
         '-m',
         'pairsamtools',
         'sort',
         '--input',
         mock_pairsam_path_2,
         '--output',
         mock_sorted_pairsam_path_2
         ],
        )

def teardown_func():
    tmpdir.cleanup()

@with_setup(setup_func, teardown_func)
def test_mock_pairsam():
    result = subprocess.check_output(
        ['python',
         '-m',
         'pairsamtools',
         'merge',
         mock_sorted_pairsam_path_1,
         mock_sorted_pairsam_path_2
         ],
        ).decode('ascii')

    # check if the two headers got transferred correctly (PG's do get modified)
    pairsam_header_1 = [l.strip() for l in open(mock_sorted_pairsam_path_1, 'r') if l.startswith('#')]
    pairsam_header_2 = [l.strip() for l in open(mock_sorted_pairsam_path_2, 'r') if l.startswith('#')]
    output_header  = [l.strip() for l in result.split('\n') if l.startswith('#')] 

    for l in pairsam_header_1:
        if not l.startswith('#@PG'):
            assert any([l in l2 for l2 in pairsam_header_1])

    for l in pairsam_header_2:
        if not l.startswith('#@PG'):
            assert any([l in l2 for l2 in pairsam_header_1])
        
    # check that all pairsam entries survived sorting:
    pairsam_body_1 = [l.strip() for l in open(mock_pairsam_path_1, 'r') 
                      if not l.startswith('#') and l.strip()]
    pairsam_body_2 = [l.strip() for l in open(mock_pairsam_path_2, 'r') 
                      if not l.startswith('#') and l.strip()]
    output_body  = [l.strip() for l in result.split('\n')
                    if not l.startswith('#') and l.strip()]

    assert len(pairsam_body_1) + len(pairsam_body_2) == len(output_body)

    # check the sorting order of the output:
    prev_pair = None
    for l in output_body:
        cur_pair = l.split('\v')[1:8]
        if prev_pair is not None:
            assert (cur_pair[0] >= prev_pair[0])
            if (cur_pair[0] == prev_pair[0]):
                assert (cur_pair[1] >= prev_pair[1])
                if (cur_pair[1] == prev_pair[1]):
                    assert (cur_pair[2] >= prev_pair[2]) 
                    if (cur_pair[2] == prev_pair[2]):
                        assert (cur_pair[3] >= prev_pair[3])

        prev_pair = cur_pair


