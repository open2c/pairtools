
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import sys
import click

import numpy as np

from collections import OrderedDict

from . import _fileio, _pairsam_format, cli, _headerops, common_io_options

UTIL_NAME = 'pairsam_stats'

@cli.command()

@click.argument(
    'input_path', 
    type=str,
    nargs=-1,
    required=False)

@click.option(
    '-o', "--output", 
    type=str, 
    default="", 
    help='output stats tsv file.')

@click.option(
    "--merge", 
    is_flag=True,
    help='If specified, merge multiple input stats files instead of calculating'
        ' statistics of a pairsam file. Merging is performed via summation of'
        ' all overlapping statistics. Non-overlapping statistics are appended to'
        ' the end of the file.',
    )

@common_io_options

def stats(input_path, output, merge, **kwargs):
    '''Calculate various statistics of a pairs/pairsam file. 

    INPUT_PATH : by default, a .pairsam file to calculate statistics.
    If not provided, the input is read from stdin.
    If --merge is specified, then INPUT_PATH is interpreted as an arbitrary number 
    of stats files to merge.
    
    The files with paths ending with .gz are gzip-decompressed. 
    '''
    stats_py(input_path, output, merge, **kwargs)

def stats_py(input_path, output, merge, **kwargs):
    if merge:
        do_merge(output, input_path, **kwargs)
        return

    instream = (_fileio.auto_open(input_path[0], mode='r', 
                                  nproc=kwargs.get('nproc_in'),
                                  command=kwargs.get('cmd_in', None)) 
                if input_path else sys.stdin)
    outstream = (_fileio.auto_open(output, mode='w', 
                                   nproc=kwargs.get('nproc_out'),
                                   command=kwargs.get('cmd_out', None)) 
                 if output else sys.stdout)


    header, body_stream = _headerops.get_header(instream)
    
    stats = OrderedDict()
    stats['total'] = 0
    stats['total_unmapped'] = 0
    stats['total_single_sided_mapped'] = 1
    stats['total_mapped'] = 0
    stats['cis'] = 0
    stats['trans'] = 0
    stats['pair_types'] = {}

    stats['cis_1kb+'] = 0
    stats['cis_2kb+'] = 0
    stats['cis_10kb+'] = 0
    stats['cis_20kb+'] = 0

    stats['chrom_freq'] = OrderedDict()
    min_log10_dist = 0
    max_log10_dist = 9
    bin_log10_spacing = 0.25
    dist_bins = (np.r_[0,
        np.round(10**np.arange(min_log10_dist, max_log10_dist+0.001, bin_log10_spacing))
        .astype(np.int)]
    )

    stats['dist_freq'] = OrderedDict([
        ('+-', np.zeros(len(dist_bins), dtype=np.int)),
        ('-+', np.zeros(len(dist_bins), dtype=np.int)),
        ('--', np.zeros(len(dist_bins), dtype=np.int)),
        ('++', np.zeros(len(dist_bins), dtype=np.int)),
        ])

    # Collecting statistics
    for line in body_stream:
        cols = line.rstrip().split(_pairsam_format.PAIRSAM_SEP)
        chrom1, pos1, strand1 = (
                cols[_pairsam_format.COL_C1], 
                int(cols[_pairsam_format.COL_P1]),
                cols[_pairsam_format.COL_S1])
        chrom2, pos2, strand2 = (
                cols[_pairsam_format.COL_C2], 
                int(cols[_pairsam_format.COL_P2]),
                cols[_pairsam_format.COL_S2])

        pair_type = cols[_pairsam_format.COL_PTYPE]
        stats['total'] += 1
        stats['pair_types'][pair_type] = stats['pair_types'].get(pair_type,0) +1
        if chrom1 == '!' and chrom2 == '!':
            stats['total_unmapped'] += 1
        elif chrom1 != '!' and chrom2 != '!':
            stats['chrom_freq'][(chrom1, chrom2)] = (
                stats['chrom_freq'].get((chrom1, chrom2),0) + 1)
            stats['total_mapped'] += 1

            if chrom1 == chrom2:
                stats['cis'] += 1
                dist = np.abs(pos2-pos1)
                bin_idx = np.searchsorted(dist_bins, dist, 'right') -1
                stats['dist_freq'][strand1+strand2][bin_idx] += 1
                if dist >= 1000:
                    stats['cis_1kb+'] += 1
                if dist >= 2000:
                    stats['cis_2kb+'] += 1
                if dist >= 10000:
                    stats['cis_10kb+'] += 1
                if dist >= 20000:
                    stats['cis_20kb+'] += 1

            else:
                stats['trans'] += 1
        else:
            stats['total_single_sided_mapped'] += 1


    # Storing statistics
    for k,v in stats.items():
        if isinstance(v, int):
            outstream.write('{}\t{}\n'.format(k,v))
        else:
            if k == 'dist_freq':
                for i in range(len(dist_bins)):
                    for dirs, freqs in v.items():
                        if i != len(dist_bins) - 1:
                            outstream.write('{}/{}-{}/{}\t{}\n'.format(
                                k, dist_bins[i], dist_bins[i+1], dirs, freqs[i])
                                )
                        else:
                            outstream.write('{}/{}+/{}\t{}\n'.format(
                                k, dist_bins[i], dirs, freqs[i]))

            if k == 'pair_types':
                for pair_type, freq in v.items():
                    outstream.write('{}/{}\t{}\n'.format(
                        k, pair_type, freq))

            if k == 'chrom_freq':
                for (chrom1, chrom2), freq in v.items():
                    outstream.write('{}/{}/{}\t{}\n'.format(
                        k, chrom1, chrom2, freq))


    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()


def do_merge(output, files_to_merge, **kwargs):

    # Parse all stats files.
    stats = []
    for stat_file in files_to_merge:
        f = _fileio.auto_open(stat_file, mode='r', 
                              nproc=kwargs.get('nproc_in'),
                              command=kwargs.get('cmd_in', None)) 
        stat = OrderedDict()
        for l in f:
            fields = l.strip().split('\t') 
            if len(fields) == 0:
                continue
            if len(fields) != 2:
                raise _fileio.ParseError(
                    '{} is not a valid stats file'.format(stat_file))
            stat[fields[0]] = int(fields[1])

        stats.append(stat)
        f.close()

    # Find a set of all possible keys. First, print overlapping keys, 
    # preserving the order. Then, add unique keys from each of the stats tables.
    out_keys = [
        k for k in stats[0]
        if all(k in stat for stat in stats)]

    for stat in stats:
        out_keys += [
            k for k in stat
            if k not in out_keys]

    # Sum all stats.
    out_stats = OrderedDict()
    for k in out_keys:
        out_stats[k] = sum(stat.get(k, 0) for stat in stats)


    # Save merged stats.
    outstream = _fileio.auto_open(output, mode='w', 
                                  nproc=kwargs.get('nproc_out'),
                                  command=kwargs.get('cmd_out', None))
    for k,v in out_stats.items():
        outstream.write('{}\t{}\n'.format(k,v))
        
    if outstream != sys.stdout:
        outstream.close()

if __name__ == '__main__':
    stats()
