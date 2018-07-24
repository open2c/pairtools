
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import sys
import click

import numpy as np

from collections import OrderedDict, Mapping

from . import _fileio, _pairsam_format, cli, _headerops, common_io_options

UTIL_NAME = 'pairtools_stats'

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
        ' statistics of a .pairs/.pairsam file. Merging is performed via summation of'
        ' all overlapping statistics. Non-overlapping statistics are appended to'
        ' the end of the file.',
    )

@common_io_options

def stats(input_path, output, merge, **kwargs):
    '''Calculate pairs statistics. 

    INPUT_PATH : by default, a .pairs/.pairsam file to calculate statistics.
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
    stats = PairCounter()

    # Collecting statistics
    for line in body_stream:
        cols = line.rstrip().split(_pairsam_format.PAIRSAM_SEP)

        # algn1:
        chrom1 = cols[_pairsam_format.COL_C1]
        pos1 = int(cols[_pairsam_format.COL_P1])
        strand1 = cols[_pairsam_format.COL_S1]
        # algn2:
        chrom2 = cols[_pairsam_format.COL_C2]
        pos2 = int(cols[_pairsam_format.COL_P2])
        strand2 = cols[_pairsam_format.COL_S2]
        # pair type:
        pair_type = cols[_pairsam_format.COL_PTYPE]

        stats.add_pair(chrom1, pos1, strand1,
                       chrom2, pos2, strand2,
                       pair_type)


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
        # use a factory method to instanciate PairCounter
        stat = PairCounter.from_file(f)
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


class PairCounter(Mapping):
    """
    A Counter for Hi-C pairs that accumulates various statistics.

    PairCounter implements two interfaces to access multi-level statistics:
    1. as a nested dict, e.g. pairCounter['pair_types']['LL']
    2. as a flat dict, with the level keys separated by '/', e.g. pairCounter['pair_types/LL']

    Other features:
    -- PairCounters can be saved into/loaded from a file
    -- multiple PairCounters can be merged via addition.
    """

    _SEP = '\t'
    _KEY_SEP = '/'

    def __init__(self, min_log10_dist=0, max_log10_dist=9, log10_dist_bin_step=0.25):
        self._stat = OrderedDict()
        # some variables used for initialization:
        # genomic distance bining for the ++/--/-+/+- distribution
        self._dist_bins = (np.r_[0,
            np.round(10**np.arange(min_log10_dist, max_log10_dist+0.001, 
                                   log10_dist_bin_step))
            .astype(np.int)]
        )

        # establish structure of an empty _stat:
        self._stat['total'] = 0
        self._stat['total_unmapped'] = 0
        self._stat['total_single_sided_mapped'] = 0
        # total_mapped = total_dups + total_nodups
        self._stat['total_mapped'] = 0
        self._stat['total_dups'] = 0
        self._stat['total_nodups'] = 0
        ########################################
        # the rest of stats are based on nodups:
        ########################################
        self._stat['cis'] = 0
        self._stat['trans'] = 0
        self._stat['pair_types'] = {}
        # to be removed:
        self._stat['dedup'] = {}

        self._stat['cis_1kb+'] = 0
        self._stat['cis_2kb+'] = 0
        self._stat['cis_4kb+'] = 0
        self._stat['cis_10kb+'] = 0
        self._stat['cis_20kb+'] = 0
        self._stat['cis_40kb+'] = 0

        self._stat['chrom_freq'] = OrderedDict()

        self._stat['dist_freq'] = OrderedDict([
            ('+-', np.zeros(len(self._dist_bins), dtype=np.int)),
            ('-+', np.zeros(len(self._dist_bins), dtype=np.int)),
            ('--', np.zeros(len(self._dist_bins), dtype=np.int)),
            ('++', np.zeros(len(self._dist_bins), dtype=np.int)),
            ])


    def __getitem__(self, key):
        if isinstance(key, str):
            # let's strip any unintentional '/'
            # from either side of the key
            key = key.strip('/')
            if self._KEY_SEP in key:
                # multi-key to access nested elements
                k_fields = key.split(self._KEY_SEP)
            else:
                # single-key access flat part of PairCounter
                # or to access highest level of hierarchy
                return self._stat[key]
        else:
            # clearly an error:
            raise ValueError(
                '{} is not a valid key: must be str'.format(key))

        # K_FIELDS:
        # process multi-key case:
        # in this case key must be in ['pair_types','chrom_freq','dist_freq','dedup']
        # get the first 'k' and keep the remainders in 'k_fields'
        k = k_fields.pop(0)
        if k in ['pair_types', 'dedup']:
            # assert there is only one element in key_fields left:
            # 'pair_types' and 'dedup' treated the same
            if len(k_fields) == 1:
                return self._stat[k][k_fields[0]]
            else:
                raise ValueError(
                    '{} is not a valid key: {} section implies 1 identifier'.format(key, k))

        elif k == 'chrom_freq':
            # assert remaining key_fields == [chr1, chr2]:
            if len(k_fields) == 2:
                return self._stat[k][tuple(k_fields)]
            else:
                raise ValueError(
                    '{} is not a valid key: {} section implies 2 identifiers'.format(key, k))

        elif k == 'dist_freq':
            # assert that last element of key_fields is the 'directions'
            # THIS IS DONE FOR CONSISTENCY WITH .stats FILE
            # SHOULD THAT BE CHANGED IN .stats AND HERE AS WELL?
            if len(k_fields) == 2:
                # assert 'dirs' in ['++','--','+-','-+']
                dirs = k_fields.pop()
                # there is only genomic distance range of the bin that's left:
                bin_range, = k_fields
                # extract left border of the bin "1000000+" or "1500-6000":
                dist_bin_left = bin_range.strip('+') if bin_range.endswith('+') \
                            else bin_range.split('-')[0]
                # get the index of that bin:
                bin_idx = np.searchsorted(self._dist_bins, int(dist_bin_left), 'right') - 1
                # store corresponding value:
                return self._stat['dist_freq'][dirs][bin_idx]
            else:
                raise ValueError(
                    '{} is not a valid key: {} section implies 2 identifiers'.format(key,k))
        else:
            raise ValueError(
                '{} is not a valid key'.format(k))


    def __iter__(self):
        return iter(self._stat)


    def __len__(self):
        return len(self._stat)


    @classmethod
    def from_file(cls, file_handle):
        """create instance of PairCounter from file

        Parameters
        ----------
        file_handle: file handle

        Returns
        -------
        PairCounter
            new PairCounter filled with the contents of the input file
        """
        # fill in from file - file_handle:
        stat_from_file = cls()
        for l in file_handle:
            fields = l.strip().split(cls._SEP) 
            if len(fields) == 0:
                # skip empty lines:
                continue
            if len(fields) != 2:
                # expect two _SEP separated values per line:
                raise _fileio.ParseError(
                    '{} is not a valid stats file'.format(file_handle.name))
            # extract key and value, then split the key:
            putative_key, putative_val =  fields[0], fields[1]
            key_fields = putative_key.split(cls._KEY_SEP)
            # we should impose a rigid structure of .stats or redo it:
            if len(key_fields)==1:
                key = key_fields[0]
                if key in stat_from_file._stat:
                    stat_from_file._stat[key] = int(fields[1])
                else:
                    raise _fileio.ParseError(
                        '{} is not a valid stats file: unknown field {} detected'.format(file_handle.name,key))
            else:
                # in this case key must be in ['pair_types','chrom_freq','dist_freq','dedup']
                # get the first 'key' and keep the remainders in 'key_fields'
                key = key_fields.pop(0)
                if key in ['pair_types', 'dedup']:
                    # assert there is only one element in key_fields left:
                    # 'pair_types' and 'dedup' treated the same
                    if len(key_fields) == 1:
                        stat_from_file._stat[key][key_fields[0]] = int(fields[1])
                    else:
                        raise _fileio.ParseError(
                            '{} is not a valid stats file: {} section implies 1 identifier'.format(file_handle.name,key))

                elif key == 'chrom_freq':
                    # assert remaining key_fields == [chr1, chr2]:
                    if len(key_fields) == 2:
                        stat_from_file._stat[key][tuple(key_fields)] = int(fields[1])
                    else:
                        raise _fileio.ParseError(
                            '{} is not a valid stats file: {} section implies 2 identifiers'.format(file_handle.name,key))

                elif key == 'dist_freq':
                    # assert that last element of key_fields is the 'directions'
                    if len(key_fields) == 2:
                        # assert 'dirs' in ['++','--','+-','-+']
                        dirs = key_fields.pop()
                        # there is only genomic distance range of the bin that's left:
                        bin_range, = key_fields
                        # extract left border of the bin "1000000+" or "1500-6000":
                        dist_bin_left = (bin_range.strip('+')
                            if bin_range.endswith('+') 
                            else bin_range.split('-')[0])
                        # get the index of that bin:
                        bin_idx = np.searchsorted(
                            stat_from_file._dist_bins, 
                            int(dist_bin_left), 'right') - 1
                        # store corresponding value:
                        stat_from_file._stat[key][dirs][bin_idx] = int(fields[1])
                    else:
                        raise _fileio.ParseError(
                            '{} is not a valid stats file: {} section implies 2 identifiers'.format(file_handle.name,key))
                else:
                    raise _fileio.ParseError(
                        '{} is not a valid stats file: unknown field {} detected'.format(file_handle.name,key))
        # return PairCounter from a non-empty dict:
        return stat_from_file


    def add_pair(self, chrom1, pos1, strand1, chrom2, pos2, strand2, pair_type):
        """Gather statistics for a Hi-C pair and add to the PairCounter.

        Parameters
        ----------
        chrom1: str
            chromosome of the first read
        pos1: int
            position of the first read
        strand1: str
            strand of the first read
        chrom2: str
            chromosome of the first read
        pos2: int
            position of the first read
        strand2: str
            strand of the first read
        pair_type: str
            type of the mapped pair of reads
        """

        self._stat['total'] += 1
        # collect pair type stats including DD:
        self._stat['pair_types'][pair_type] = self._stat['pair_types'].get(pair_type,0) + 1
        if chrom1 == '!' and chrom2 == '!':
            self._stat['total_unmapped'] += 1
        elif chrom1 != '!' and chrom2 != '!':
            self._stat['total_mapped'] += 1
            # only mapped ones can be duplicates:
            if pair_type == 'DD':
                self._stat['total_dups'] += 1
            else:
                self._stat['total_nodups'] += 1
                self._stat['chrom_freq'][(chrom1, chrom2)] = (
                    self._stat['chrom_freq'].get((chrom1, chrom2), 0) + 1)

                if chrom1 == chrom2:
                    self._stat['cis'] += 1
                    dist = np.abs(pos2-pos1)
                    bin_idx = np.searchsorted(self._dist_bins, dist, 'right') - 1
                    self._stat['dist_freq'][strand1+strand2][bin_idx] += 1
                    if dist >= 1000:
                        self._stat['cis_1kb+'] += 1
                    if dist >= 2000:
                        self._stat['cis_2kb+'] += 1
                    if dist >= 4000:
                        self._stat['cis_4kb+'] += 1
                    if dist >= 10000:
                        self._stat['cis_10kb+'] += 1
                    if dist >= 20000:
                        self._stat['cis_20kb+'] += 1
                    if dist >= 40000:
                        self._stat['cis_40kb+'] += 1

                else:
                    self._stat['trans'] += 1
        else:
            self._stat['total_single_sided_mapped'] += 1


    def __add__(self, other):
        # both PairCounter are implied to have a list of common fields:
        #
        # 'total', 'total_unmapped', 'total_single_sided_mapped', 'total_mapped',
        # 'cis', 'trans', 'pair_types', 'cis_1kb+', 'cis_2kb+',
        # 'cis_10kb+', 'cis_20kb+', 'chrom_freq', 'dist_freq', 'dedup'
        #
        # initialize empty PairCounter for the result of summation:
        sum_stat = PairCounter()
        # use the empty PairCounter to iterate over:
        for k,v in sum_stat._stat.items():
            # not nested fields are summed trivially:
            if isinstance(v, int):
                sum_stat._stat[k] = self._stat[k] + other._stat[k]
            # sum nested dicts/arrays in a context dependet manner:
            else:
                if k in ['pair_types','dedup']:
                    # handy function for summation of a pair of dicts:
                    # https://stackoverflow.com/questions/10461531/merge-and-sum-of-two-dictionaries
                    sum_dicts = lambda dict_x,dict_y: { 
                                    key: dict_x.get(key, 0) + dict_y.get(key, 0)
                                            for key in set(dict_x) | set(dict_y) }
                    # sum a pair of corresponding dicts:
                    sum_stat._stat[k] = sum_dicts(self._stat[k], other._stat[k])
                if k == 'chrom_freq':
                    # union list of keys (chr1,chr2) with potential duplicates:
                    union_keys_with_dups = list(self._stat[k].keys()) + list(other._stat[k].keys())
                    # OrderedDict.fromkeys will take care of keys' order and duplicates in a consistent manner:
                    # https://stackoverflow.com/questions/1720421/how-to-concatenate-two-lists-in-python
                    # last comment to the 3rd Answer
                    sum_stat._stat[k] = OrderedDict.fromkeys( union_keys_with_dups )
                    # perform a summation:
                    for union_key in sum_stat._stat[k]:
                        sum_stat._stat[k][union_key] = (
                            self._stat[k].get(union_key, 0) 
                            + other._stat[k].get(union_key, 0)
                        )
                if k == 'dist_freq':
                    for dirs in sum_stat[k]:
                        sum_stat[k][dirs] = self._stat[k][dirs] + other._stat[k][dirs]
        return sum_stat


    # we need this to be able to sum(list_of_PairCounters)
    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)


    def flatten(self):
        """return a flattened OrderedDic (formatted same way as .stats file)

        """
        # OrderedDict for flat store:
        flat_stat = OrderedDict()

        # Storing statistics
        for k,v in self._stat.items():
            if isinstance(v, int):
                flat_stat[k] = v
            # store nested dicts/arrays in a context dependet manner:
            # nested categories are stored only if they are non-trivial
            else:
                if (k == 'dist_freq') and v:
                    for i in range(len(self._dist_bins)):
                        for dirs, freqs in v.items():
                            # last bin is treated differently: "100000+" vs "1200-3000":
                            if i != len(self._dist_bins) - 1:
                                formatted_key = self._KEY_SEP.join(['{}','{}-{}','{}']).format(
                                                k, self._dist_bins[i], self._dist_bins[i+1], dirs)
                            else:
                                formatted_key = self._KEY_SEP.join(['{}','{}+','{}']).format(
                                                k, self._dist_bins[i], dirs)
                            #store key,value pair:
                            flat_stat[formatted_key] = freqs[i]
                elif (k in ['pair_types','dedup']) and v:
                    # 'pair_types' and 'dedup' are simple dicts inside,
                    # treat them the exact same way:
                    for k_item, freq in v.items():
                        formatted_key = self._KEY_SEP.join(['{}','{}']).format(k, k_item)
                        #store key,value pair:
                        flat_stat[formatted_key] = freq
                elif (k == 'chrom_freq') and v:
                    for (chrom1, chrom2), freq in v.items():
                        formatted_key = self._KEY_SEP.join(['{}','{}','{}']).format(k, chrom1, chrom2)
                        #store key,value pair:
                        flat_stat[formatted_key] = freq

        # return flattened OrderedDict
        return flat_stat


    def save(self, outstream):
        """save PairCounter to tab-delimited text file.
        Flattened version of PairCounter is stored in the file.

        Parameters
        ----------
        outstream: file handle


        Note
        ----
        The order of the keys is not guaranteed
        Merging several .stats is not associative with respect to key order:
        merge(A,merge(B,C)) != merge(merge(A,B),C).

        Theys shou5ld match exactly, however, when soprted:
        sort(merge(A,merge(B,C))) == sort(merge(merge(A,B),C))
        """

        # write flattened version of the PairCounter to outstream
        for k,v in self.flatten().items():
            outstream.write('{}{}{}\n'.format(k, self._SEP, v))


if __name__ == '__main__':
    stats()

