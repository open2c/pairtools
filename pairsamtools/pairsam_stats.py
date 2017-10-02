
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
    
    The files with paths ending with .gz/.lz4 are decompressed by pbgzip/lz4c. 
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


    # new stats class stuff would come here ...
    stats = StatObject()

    # Collecting statistics
    for line in body_stream:
        cols = line.rstrip().split(_pairsam_format.PAIRSAM_SEP)

        algn1 = {
            'chrom': cols[_pairsam_format.COL_C1],
            'pos': int(cols[_pairsam_format.COL_P1]),
            'strand': cols[_pairsam_format.COL_S1] }

        algn2 = {
            'chrom': cols[_pairsam_format.COL_C2],
            'pos': int(cols[_pairsam_format.COL_P2]),
            'strand': cols[_pairsam_format.COL_S2] }

        pair_type = cols[_pairsam_format.COL_PTYPE]
        #
        #
        stats.update(algn1, algn2, pair_type)



    # save statistics to file ...
    stats.save(outstream)



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
        stat = StatObject(file_handle = f)
        stats.append(stat)
        f.close()

    # combine stats from several files (files_to_merge):
    out_stat = sum(stats)


    # Save merged stats.
    outstream = _fileio.auto_open(output, mode='w', 
                                  nproc=kwargs.get('nproc_out'),
                                  command=kwargs.get('cmd_out', None))

    # save statistics to file ...
    out_stat.save(outstream)
        
    if outstream != sys.stdout:
        outstream.close()


class StatObject(object):
    """docstring for StatObject:
    A class that defines a simple container
    for sequencing statistics and simple operations
    for the class (accumulation, merge and output)."""
    def __init__(self, file_handle = None, stat_dict = None):

        if stat_dict:
            self.stat = stat_dict
        else:
            # initialize empty stat dictironary:
            self.stat = OrderedDict()
            # fill it from file or establish an empty one:
            if file_handle:
                # fill in from file_handle:
                for l in file_handle:
                    fields = l.strip().split('\t') 
                    if len(fields) == 0:
                        continue
                    if len(fields) != 2:
                        raise _fileio.ParseError(
                            '{} is not a valid stats file'.format(file_handle.name))
                    self.stat[fields[0]] = int(fields[1])
            else:
                # some more variables used for initialization:
                self.min_log10_dist = 0
                self.max_log10_dist = 9
                self.bin_log10_spacing = 0.25
                self.dist_bins = (np.r_[0,
                    np.round(10**np.arange(self.min_log10_dist, self.max_log10_dist+0.001, self.bin_log10_spacing))
                    .astype(np.int)]
                )

                # establish structure of an empty stat:
                self.stat['total'] = 0
                self.stat['total_unmapped'] = 0
                self.stat['total_single_sided_mapped'] = 1
                self.stat['total_mapped'] = 0
                self.stat['cis'] = 0
                self.stat['trans'] = 0
                self.stat['pair_types'] = {}

                self.stat['cis_1kb+'] = 0
                self.stat['cis_2kb+'] = 0
                self.stat['cis_10kb+'] = 0
                self.stat['cis_20kb+'] = 0

                self.stat['chrom_freq'] = OrderedDict()

                self.stat['dist_freq'] = OrderedDict([
                    ('+-', np.zeros(len(self.dist_bins), dtype=np.int)),
                    ('-+', np.zeros(len(self.dist_bins), dtype=np.int)),
                    ('--', np.zeros(len(self.dist_bins), dtype=np.int)),
                    ('++', np.zeros(len(self.dist_bins), dtype=np.int)),
                    ])


    def update(self, algn1, algn2, pair_type):
        # extract chrom, position and strand from each of the alignmentns:
        chrom1, pos1, strand1 = ( algn1['chrom'], algn1['pos'], algn1['strand'] )
        chrom2, pos2, strand2 = ( algn2['chrom'], algn2['pos'], algn2['strand'] )

        self.stat['total'] += 1
        self.stat['pair_types'][pair_type] = self.stat['pair_types'].get(pair_type,0) + 1
        if chrom1 == '!' and chrom2 == '!':
            self.stat['total_unmapped'] += 1
        elif chrom1 != '!' and chrom2 != '!':
            self.stat['chrom_freq'][(chrom1, chrom2)] = (
                self.stat['chrom_freq'].get((chrom1, chrom2),0) + 1)
            self.stat['total_mapped'] += 1

            if chrom1 == chrom2:
                self.stat['cis'] += 1
                dist = np.abs(pos2-pos1)
                bin_idx = np.searchsorted(self.dist_bins, dist, 'right') - 1
                self.stat['dist_freq'][strand1+strand2][bin_idx] += 1
                if dist >= 1000:
                    self.stat['cis_1kb+'] += 1
                if dist >= 2000:
                    self.stat['cis_2kb+'] += 1
                if dist >= 10000:
                    self.stat['cis_10kb+'] += 1
                if dist >= 20000:
                    self.stat['cis_20kb+'] += 1

            else:
                self.stat['trans'] += 1
        else:
            self.stat['total_single_sided_mapped'] += 1


    def __add__(self, other):
        stats = [self.stat, other.stat]

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

        # return the result of a merger:
        return StatObject(stat_dict = out_stats)


    # we need this to be able to sum(list_of_StatObjects)
    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)


    def save(self, outstream):
        # Storing statistics
        for k,v in self.stat.items():
            # this should work for the stat initialized with the file:
            if isinstance(v, int):
                outstream.write('{}\t{}\n'.format(k,v))
            # this is used to store newly initialized stat-object:
            else:
                if k == 'dist_freq':
                    for i in range(len(self.dist_bins)):
                        for dirs, freqs in v.items():
                            if i != len(self.dist_bins) - 1:
                                outstream.write('{}/{}-{}/{}\t{}\n'.format(
                                    k, self.dist_bins[i], self.dist_bins[i+1], dirs, freqs[i])
                                    )
                            else:
                                outstream.write('{}/{}+/{}\t{}\n'.format(
                                    k, self.dist_bins[i], dirs, freqs[i]))

                if k == 'pair_types':
                    for pair_type, freq in v.items():
                        outstream.write('{}/{}\t{}\n'.format(
                            k, pair_type, freq))

                if k == 'chrom_freq':
                    for (chrom1, chrom2), freq in v.items():
                        outstream.write('{}/{}/{}\t{}\n'.format(
                            k, chrom1, chrom2, freq))


if __name__ == '__main__':
    stats()









