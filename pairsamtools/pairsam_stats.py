
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import sys
import click

import numpy as np

from collections import OrderedDict, Mapping

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
        stat = StatObject.from_file(f)
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


class StatObject(Mapping):
    """
    Container for various sequencing statistics

    New object can be created as empty, from file and from existing dict.
    Multiple instances of StatsObject can be merged using summation.
    Instance of StatsObject can be saved to a text file.


    Parameters
    ----------
    stat_dict: dict-like, optional
        dictionary used to initialize StatObject

    Methods
    -------
    from_file(file_handle)
        create StatsObject from file

        file_handle: file handle

    update(algn1, algn2, pair_type)
        update existing StatObject object with
        information from a mapped paired-end read

        algn1: tuple-like
            tuple contents: (chrom, pos, strand)
        algn2: tuple-like
            tuple contents: (chrom, pos, strand)
        pair_type: str
            type of the mapped pair

        Note
        ----
        this method can update StatObject
        created as StatObject() only.
        StatObject created using from_file(),
        or using StatObject(dict) cannot be updated now.

    save(outstream)
        save to .stats text file

        outstream: file handle


    """
    def __init__(self, stat_dict = None):
        if stat_dict:
            # init from non-empty dictionary
            # should we deep_copy instead ?
            self._stat = stat_dict
        else:
            # create empty dict
            self._stat = OrderedDict()
            # some variables used for initialization:
            # genomic distnace bining for the ++/--/-+/+- distribution
            self._min_log10_dist = 0
            self._max_log10_dist = 9
            self._bin_log10_spacing = 0.25
            self._dist_bins = (np.r_[0,
                np.round(10**np.arange(self._min_log10_dist, self._max_log10_dist+0.001, self._bin_log10_spacing))
                .astype(np.int)]
            )

            # establish structure of an empty _stat:
            self._stat['total'] = 0
            self._stat['total_unmapped'] = 0
            self._stat['total_single_sided_mapped'] = 0
            self._stat['total_mapped'] = 0
            self._stat['cis'] = 0
            self._stat['trans'] = 0
            self._stat['pair_types'] = {}

            self._stat['cis_1kb+'] = 0
            self._stat['cis_2kb+'] = 0
            self._stat['cis_10kb+'] = 0
            self._stat['cis_20kb+'] = 0

            self._stat['chrom_freq'] = OrderedDict()

            self._stat['dist_freq'] = OrderedDict([
                ('+-', np.zeros(len(self._dist_bins), dtype=np.int)),
                ('-+', np.zeros(len(self._dist_bins), dtype=np.int)),
                ('--', np.zeros(len(self._dist_bins), dtype=np.int)),
                ('++', np.zeros(len(self._dist_bins), dtype=np.int)),
                ])


    def __getitem__(self, key):
        return self._stat[key]

    def __iter__(self):
        return iter(self._stat)

    def __len__(self):
        return len(self._stat)


    def from_file(self, file_handle):
        """create instance of StatObject from file

        Parameters
        ----------
        file_handle: file handle

        Returns
        -------
        StatObject
            new instance of StatObject
            filled with the contents of
            the input file

        Note
        ----
        instance of StatObject returned
        by this method cannot be updated 
        since it is a flat version of a dict
        """
        # fill in from file - file_handle:
        stat_from_file = OrderedDict()
        for l in file_handle:
            fields = l.strip().split('\t') 
            if len(fields) == 0:
                continue
            if len(fields) != 2:
                raise _fileio.ParseError(
                    '{} is not a valid stats file'.format(file_handle.name))
            stat_from_file[fields[0]] = int(fields[1])
        # return StatObject from a non-empty dict:
        return StatObject(stat_from_file)



    def update(self, algn1, algn2, pair_type):
        """update existing StatObject with info from mapped read

        Parameters
        ----------
        algn1: tuple-like
            tuple contents: (chrom, pos, strand)
        algn2: tuple-like
            tuple contents: (chrom, pos, strand)
        pair_type: str
            type of the mapped pair
            e.g. CX,LL,MN,NN, etc.


        Note
        ----
        Only instances of StatObject
        created as StatObject() (i.e., as empty)
        can be updated by this method
        It all has to do with some dicts being flat,
        while others nested.
        This will be addressed in the future.
        """

        # extract chrom, position and strand from each of the alignmentns:
        chrom1, pos1, strand1 = ( algn1['chrom'], algn1['pos'], algn1['strand'] )
        chrom2, pos2, strand2 = ( algn2['chrom'], algn2['pos'], algn2['strand'] )

        self._stat['total'] += 1
        self._stat['pair_types'][pair_type] = self._stat['pair_types'].get(pair_type,0) + 1
        if chrom1 == '!' and chrom2 == '!':
            self._stat['total_unmapped'] += 1
        elif chrom1 != '!' and chrom2 != '!':
            self._stat['chrom_freq'][(chrom1, chrom2)] = (
                self._stat['chrom_freq'].get((chrom1, chrom2),0) + 1)
            self._stat['total_mapped'] += 1

            if chrom1 == chrom2:
                self._stat['cis'] += 1
                dist = np.abs(pos2-pos1)
                bin_idx = np.searchsorted(self._dist_bins, dist, 'right') - 1
                self._stat['dist_freq'][strand1+strand2][bin_idx] += 1
                if dist >= 1000:
                    self._stat['cis_1kb+'] += 1
                if dist >= 2000:
                    self._stat['cis_2kb+'] += 1
                if dist >= 10000:
                    self._stat['cis_10kb+'] += 1
                if dist >= 20000:
                    self._stat['cis_20kb+'] += 1

            else:
                self._stat['trans'] += 1
        else:
            self._stat['total_single_sided_mapped'] += 1


    def __add__(self, other):
        stats = [self._stat, other._stat]

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
        return StatObject(out_stats)


    # we need this to be able to sum(list_of_StatObjects)
    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)


    def save(self, outstream):
        """save StatObject to tab-delimited text file

        Parameters
        ----------
        outstream: file handle


        Note
        ----
        The order of the keys is not guaranteed
        Merging several .stats is not associative with respect to key order:
        merge(A,merge(B,C)) != merge(merge(A,B),C).

        Theys should match exactly, however, when soprted:
        sort(merge(A,merge(B,C))) == sort(merge(merge(A,B),C))
        """

        # Storing statistics
        for k,v in self._stat.items():
            # this should work for the stat initialized with the file:
            if isinstance(v, int):
                outstream.write('{}\t{}\n'.format(k,v))
            # this is used to store newly initialized stat-object:
            else:
                if k == 'dist_freq':
                    for i in range(len(self._dist_bins)):
                        for dirs, freqs in v.items():
                            if i != len(self._dist_bins) - 1:
                                outstream.write('{}/{}-{}/{}\t{}\n'.format(
                                    k, self._dist_bins[i], self._dist_bins[i+1], dirs, freqs[i])
                                    )
                            else:
                                outstream.write('{}/{}+/{}\t{}\n'.format(
                                    k, self._dist_bins[i], dirs, freqs[i]))

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









