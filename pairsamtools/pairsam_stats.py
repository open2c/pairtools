
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import sys
import click

import numpy as np

from collections import OrderedDict

from . import _common, cli, _headerops

UTIL_NAME = 'pairsam_stats'

@cli.command()

@click.argument(
    'pairsam_path', 
    type=str,
    required=False)

@click.option(
    '-o', "--output", 
    type=str, 
    default="", 
    help='output stats tsv file.')

@click.option(
    '-m', "--merge", 
    type=str, 
    default="", 
    multiple=True,
    help='If provided, merge multiple stats files instead of calculating'
        ' statistics on the input file. Provide each of the files to merge using'
        ' a separate -m flag.',
    )


def stats(pairsam_path, output, merge):
    '''Calculate various statistics of a pairs/pairsam file. 

    PAIRSAM_PATH : input .pairsam file. If the path ends with .gz, the input is
    gzip-decompressed. By default, the input is read from stdin.
    '''

    if merge:
        do_merge(output, merge)
        return

    instream = (_common.open_bgzip(pairsam_path, mode='r') 
                if pairsam_path else sys.stdin)
    outstream = (_common.open_bgzip(output, mode='w') 
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
    stats['chrom_freq'] = OrderedDict()
    min_log10_dist = 0
    max_log10_dist = 9
    bin_log10_spacing = 0.25
    dist_bins = (np.r_[0,
        np.round(10**np.arange(min_log10_dist, max_log10_dist+0.001, bin_log10_spacing))
        .astype(np.int)]
    )

    stats['dist_bins'] = dist_bins
    stats['dist_freq'] = OrderedDict([
        ('+-', np.zeros(len(dist_bins), dtype=np.int)),
        ('-+', np.zeros(len(dist_bins), dtype=np.int)),
        ('--', np.zeros(len(dist_bins), dtype=np.int)),
        ('++', np.zeros(len(dist_bins), dtype=np.int)),
        ])

    for line in body_stream:
        cols = line[:-1].split(_common.PAIRSAM_SEP)
        chrom1, pos1, strand1 = (
                cols[_common.COL_C1], 
                int(cols[_common.COL_P1]),
                cols[_common.COL_S1])
        chrom2, pos2, strand2 = (
                cols[_common.COL_C2], 
                int(cols[_common.COL_P2]),
                cols[_common.COL_S2])

        pair_type = cols[_common.COL_PTYPE]
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
                bin_idx = np.searchsorted(dist_bins, pos2-pos1, 'right') -1
                stats['dist_freq'][strand1+strand2][bin_idx] += 1

            else:
                stats['trans'] += 1
        else:
            stats['total_single_sided_mapped'] += 1

    for k,v in stats.items():
        if isinstance(v, int):
            outstream.write('{}\t{}\n'.format(k,v))
        else:
            if k == 'dist_freq':
                for dirs, freqs in v.items():
                    for i, freq in enumerate(freqs):
                        if i != len(freqs) -1:
                            outstream.write('{}/{}/{}-{}\t{}\n'.format(
                                k, dirs, dist_bins[i], dist_bins[i+1], freq)
                                )
                        else:
                            outstream.write('{}/{}/{}+\t{}\n'.format(
                                k, dirs, dist_bins[i], freq))

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


def do_merge(output, files_to_merge):

    # Parse all stats files.
    stats = []
    for stat_file in files_to_merge:
        f = _common.open_bgzip(stat_file,'r')
        stat = OrderedDict()
        for l in f:
            fields = l.strip().split('\t') 
            if len(fields) == 0:
                continue
            if len(fields) != 2:
                raise Exception(
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
    outstream = (_common.open_bgzip(output, mode='w') 
                 if output else sys.stdout)

    for k,v in out_stats.items():
        outstream.write('{}\t{}\n'.format(k,v))
        
    if outstream != sys.stdout:
        outstream.close()

if __name__ == '__main__':
    stats()
