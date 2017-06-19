# -*- coding: utf-8 -*-
import os
import sys
import subprocess
from nose.tools import assert_raises

testdir = os.path.dirname(os.path.realpath(__file__))
mock_pairsam_path = os.path.join(testdir, 'data', 'mock.pairsam')


def test_preserve():
    try:
        result = subprocess.check_output(
            ['python',
             '-m',
             'pairtools',
             'select',
             'True',
             mock_pairsam_path],
            ).decode('ascii')
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

        
    pairsam_body = [l.strip() for l in open(mock_pairsam_path, 'r') 
                    if not l.startswith('#') and l.strip()]
    output_body  = [l.strip() for l in result.split('\n')
                    if not l.startswith('#') and l.strip()]
    assert all(l in pairsam_body for l in output_body)


def test_equal():
    try:
        result = subprocess.check_output(
            ['python',
             '-m',
             'pairtools',
             'select',
             '(pair_type == "CX") or (pair_type == "LL")',
             mock_pairsam_path],
            ).decode('ascii')
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e
    print(result)

    pairsam_body = [l.strip() for l in open(mock_pairsam_path, 'r') 
                    if not l.startswith('#') and l.strip()]
    output_body  = [l.strip() for l in result.split('\n')
                    if not l.startswith('#') and l.strip()]

    assert all(l.split('\t')[7] in ['CX', 'LL'] for l in output_body)
    assert all(l in output_body
               for l in pairsam_body 
               if l.split('\t')[7] in ['CX', 'LL'])


def test_csv():
    try:
        result = subprocess.check_output(
            ['python',
             '-m',
             'pairtools',
             'select',
             'csv_match(pair_type, "CX,LL")',
             mock_pairsam_path],
            ).decode('ascii')
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e
    print(result)

    pairsam_body = [l.strip() for l in open(mock_pairsam_path, 'r') 
                    if not l.startswith('#') and l.strip()]
    output_body  = [l.strip() for l in result.split('\n')
                    if not l.startswith('#') and l.strip()]

    assert all(l.split('\t')[7] in ['CX', 'LL'] for l in output_body)
    assert all(l in output_body
               for l in pairsam_body 
               if l.split('\t')[7] in ['CX', 'LL'])


def test_wildcard():
    try:
        result = subprocess.check_output(
            ['python',
             '-m',
             'pairtools',
             'select',
             'wildcard_match(pair_type, "*L")',
             mock_pairsam_path],
            ).decode('ascii')
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e
    print(result)

    pairsam_body = [l.strip() for l in open(mock_pairsam_path, 'r') 
                    if not l.startswith('#') and l.strip()]
    output_body  = [l.strip() for l in result.split('\n')
                    if not l.startswith('#') and l.strip()]

    assert all(l.split('\t')[7] in ['NL', 'ML', 'CL', 'LL'] for l in output_body)
    assert all(l in output_body
               for l in pairsam_body 
               if l.split('\t')[7] in ['NL', 'ML', 'CL', 'LL'])


def test_regex():
    try:
        result = subprocess.check_output(
            ['python',
             '-m',
             'pairtools',
             'select',
             'regex_match(pair_type, "[NM]L")',
             mock_pairsam_path],
            ).decode('ascii')
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e
    print(result)

    pairsam_body = [l.strip() for l in open(mock_pairsam_path, 'r') 
                    if not l.startswith('#') and l.strip()]
    output_body  = [l.strip() for l in result.split('\n')
                    if not l.startswith('#') and l.strip()]

    assert all(l.split('\t')[7] in ['NL', 'ML'] for l in output_body)
    assert all(l in output_body
               for l in pairsam_body 
               if l.split('\t')[7] in ['NL', 'ML'])
    
