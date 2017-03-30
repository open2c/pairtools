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
             'pairsamtools',
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
             'pairsamtools',
             'select',
             '(PAIR_TYPE == "CX") or (PAIR_TYPE == "LL")',
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
             'pairsamtools',
             'select',
             'csv_match(PAIR_TYPE, "CX,LL")',
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
             'pairsamtools',
             'select',
             'wildcard_match(PAIR_TYPE, "*L")',
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
             'pairsamtools',
             'select',
             'regex_match(PAIR_TYPE, "[NM]L")',
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
    
