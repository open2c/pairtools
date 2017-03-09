# -*- coding: utf-8 -*-
import sys
import numpy as np
import classify_sam

from nose.tools import assert_raises

def test_python_version():
    assert (sys.version_info[0] == 3), 'Use Python 3!'

def test_tests():
    assert (classify_sam.parse_cigar(b'50M') == (0,0,50,50))

def test_tests():
    assert (classify_sam.parse_cigar(b'50M') == (0,0,50,50))
