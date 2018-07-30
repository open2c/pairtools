# -*- coding: utf-8 -*-
import os
import sys

from nose.tools import assert_raises

import subprocess

testdir = os.path.dirname(os.path.realpath(__file__))

from pairtools import parse, parse_algn, parse_cigar

def test_python_version():
    assert (sys.version_info[0] == 3), 'Use Python 3!'

def test_parse_cigar():
    assert (parse_cigar('*') == 
        {
         'cigar' : '*',
         'read_len': 0, 
         'matched_bp': 0, 
         'algn_ref_span': 0, 
         'algn_read_span': 0, 
         'clip5_ref': 0, 
         'clip3_ref': 0})

    assert (parse_cigar('50M') == 
        {
         'cigar' : '50M',
         'read_len': 50, 
         'matched_bp': 50, 
         'algn_ref_span': 50, 
         'algn_read_span': 50, 
         'clip5_ref': 0, 
         'clip3_ref': 0})

    assert (parse_cigar('40M10S') == 
        {
         'cigar' : '40M10S',
         'read_len': 50, 
         'matched_bp': 40, 
         'algn_ref_span': 40, 
         'algn_read_span': 40, 
         'clip5_ref': 0, 
         'clip3_ref': 10})

    assert (parse_cigar('10S40M') == 
        {
         'cigar' : '10S40M',
         'read_len': 50, 
         'matched_bp': 40, 
         'algn_ref_span': 40, 
         'algn_read_span': 40, 
         'clip5_ref': 10, 
         'clip3_ref': 0})

    assert (parse_cigar('10S30M10S') == 
        {
         'cigar' : '10S30M10S',
         'read_len': 50, 
         'matched_bp': 30, 
         'algn_ref_span': 30, 
         'algn_read_span': 30, 
         'clip5_ref': 10, 
         'clip3_ref': 10})

    assert (parse_cigar('30M10I10M') == 
        {
         'cigar' : '30M10I10M',
         'read_len': 50, 
         'matched_bp': 40, 
         'algn_ref_span': 40, 
         'algn_read_span': 50, 
         'clip5_ref': 0, 
         'clip3_ref': 0})

    assert (parse_cigar('30M10D10M10S') == 
        {
         'cigar' : '30M10D10M10S',
         'read_len': 50, 
         'matched_bp': 40, 
         'algn_ref_span': 50, 
         'algn_read_span': 40, 
         'clip5_ref': 0, 
         'clip3_ref': 10})


def test_parse_algn():
    min_mapq = 50

    sam='SRR1658570.5\t65\tchr12\t24316205\t60\t90M11S\t'
    '=\t46893391\t22577187\t.\t.\t'
    'NM:i:1\tMD:Z:36A53\tAS:i:85\tXS:i:19'
    samcols = sam.split('\t')
    parsed_algn = parse_algn(samcols, min_mapq)
    assert parsed_algn == {
         'chrom': 'chr12', 
         'pos': 24316205, 
         'pos5': 24316205, 
         'pos3': 24316294, 
         'pos': 24316205, 
         'strand': '+', 
         'dist_to_5': 0, 
         'dist_to_3': 11, 
         'mapq': 60, 
         'is_unique': True, 
         'is_mapped': True, 
         'is_linear': True, 
         'cigar' : '90M11S',
         'algn_ref_span': 90, 
         'algn_read_span': 90,
         'matched_bp': 90, 
         'clip3_ref': 11, 
         'clip5_ref': 0, 
         'read_len': 101,
         'type':'U'}

    sam = ('readid01\t65\tchr1\t10\t60\t50M\tchr1\t200\t0\tSEQ\tPHRED'
          '\tFLAG1\tFLAG2\tSIMULATED:readid01,chr1,chr1,10,200,+,+,UU')
    samcols = sam.split('\t')
    parsed_algn = parse_algn(samcols, min_mapq, True)
    assert parsed_algn == {
         'chrom': 'chr1', 
         'pos': 59, 
         'pos5': 10, 
         'pos3': 59, 
         'strand': '+', 
         'dist_to_5': 0, 
         'dist_to_3': 0, 
         'mapq': 60, 
         'is_unique': True, 
         'is_mapped': True, 
         'is_linear': True, 
         'cigar' : '50M',
         'algn_ref_span': 50, 
         'algn_read_span': 50,
         'matched_bp': 50, 
         'clip3_ref': 0,
         'clip5_ref': 0, 
         'read_len': 50,
         'type':'U'}


    sam = ('readid10\t77\t*\t0\t0\t*\t*\t0\t0\tSEQ\tPHRED'
           '\tFLAG1\tFLAG2\tSIMULATED:readid10,!,!,0,0,-,-,NN')
    samcols = sam.split('\t')
    parsed_algn = parse_algn(samcols, min_mapq)
    assert parsed_algn == {
         'chrom': '!', 
         'pos': 0, 
         'pos5': 0, 
         'pos3': 0, 
         'strand': '-', 
         'dist_to_5': 0, 
         'dist_to_3': 0, 
         'mapq': 0, 
         'is_unique': False, 
         'is_mapped': False, 
         'is_linear': True, 
         'cigar' : '*',
         'algn_ref_span': 0, 
         'algn_read_span': 0,
         'matched_bp': 0, 
         'clip3_ref': 0,
         'clip5_ref': 0, 
         'read_len': 0,
         'type':'N'}


def test_mock_sam():
    mock_sam_path = os.path.join(testdir, 'data', 'mock.sam')
    mock_chroms_path = os.path.join(testdir, 'data', 'mock.chrom.sizes')
    try:
        result = subprocess.check_output(
            ['python',
             '-m',
             'pairtools',
             'parse',
             '--walks-policy',
             'mask',
             '-c',
             mock_chroms_path,
             mock_sam_path],
            ).decode('ascii')
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    # check if the header got transferred correctly
    sam_header = [l.strip() for l in open(mock_sam_path, 'r') if l.startswith('@')]
    pairsam_header = [l.strip() for l in result.split('\n') if l.startswith('#')] 
    for l in sam_header:
        assert any([l in l2 for l2 in pairsam_header])

    # check that the pairs got assigned properly
    for l in result.split('\n'):
        if l.startswith('#') or not l:
            continue

        assigned_pair = l.split('\t')[1:8]
        simulated_pair = l.split('SIMULATED:',1)[1].split('\031',1)[0].split(',')
        print(assigned_pair)
        print(simulated_pair)
        print()

        assert assigned_pair == simulated_pair

