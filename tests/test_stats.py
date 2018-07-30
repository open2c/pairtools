# -*- coding: utf-8 -*-
import os
import sys
import subprocess
from nose.tools import assert_raises

testdir = os.path.dirname(os.path.realpath(__file__))

def test_mock_pairsam():
    mock_pairsam_path = os.path.join(testdir, 'data', 'mock.pairsam')
    try:
        result = subprocess.check_output(
            ['python',
             '-m',
             'pairtools',
             'stats',
             mock_pairsam_path],
            ).decode('ascii')
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

        
    stats = dict(l.strip().split('\t') 
                 for l in result.split('\n')
                if not l.startswith('#') and l.strip())

    for k in stats:
        stats[k] = int(stats[k])
    print(stats)

    assert stats['total'] == 8
    assert stats['total_single_sided_mapped'] == 2
    assert stats['total_mapped'] == 5
    assert stats['cis'] == 3
    assert stats['trans'] == 2
    assert stats['pair_types/UU'] == 4
    assert stats['pair_types/NU'] == 1
    assert stats['pair_types/WW'] == 1
    assert stats['pair_types/UR'] == 1
    assert stats['pair_types/MU'] == 1
    assert stats['chrom_freq/chr1/chr2'] ==  1
    assert stats['chrom_freq/chr1/chr1'] ==  3
    assert stats['chrom_freq/chr2/chr3'] ==  1
    assert all(stats[k]==0 
               for k in stats
               if k.startswith('dist_freq')
               and k not in ['dist_freq/1-2/++',
                             'dist_freq/2-3/++',
                             'dist_freq/32-56/++'])

    assert stats['dist_freq/1-2/++'] == 1
    assert stats['dist_freq/2-3/++'] == 1
    assert stats['dist_freq/32-56/++'] == 1

