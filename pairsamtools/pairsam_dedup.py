#!/usr/bin/env python
# -*- coding: utf-8  -*-
import sys
import ast 
import warnings

import click

import numpy as np

from . import _dedup, _io, _pairsam_format, _headerops, cli


UTIL_NAME = 'pairsam_dedup'

# you don't need to load more than 10k lines at a time b/c you get out of the 
# CPU cache, so this parameter is not adjustable
MAX_LEN = 10000 


@cli.command()
@click.argument(
    'pairsam_path', 
    type=str,
    required=False)
@click.option(
    "-o", "--output", 
    type=str, 
    default="", 
    help='output file for pairs after duplicate removal.'
        ' If the path ends with .gz, the output is bgzip-compressed.'
        ' By default, the output is printed into stdout.')
@click.option(
    "--output-dups",
    type=str, 
    default="", 
    help='output file for duplicates. '
        ' If the path ends with .gz, the output is bgzip-compressed.'
        ' By default, duplicates are dropped.')
@click.option(
    "--stats-file", 
    type=str, 
    default="", 
    help='output file for duplicate statistics. '
        ' If file exists, it will be open in the append mode.'
        ' If the path ends with .gz, the output is bgzip-compressed.'
        ' By default, statistics are not printed.')
@click.option(
    "--max-mismatch",
    type=int, 
    default=3,
    help='Pairs with both sides mapped within this distance (bp) from each '
         'other are considered duplicates.')
@click.option(
    '--method',
    type=click.Choice(['max', 'sum']),
    default="max",  
    help='define the mismatch as either the max or the sum of the mismatches of'
        'the genomic locations of the both sides of the two compared molecules',
    show_default=True,)
@click.option(
    "--sep",
    type=str, 
    default=_pairsam_format.PAIRSAM_SEP_ESCAPE, 
    help=r"Separator (\t, \v, etc. characters are "
          "supported, pass them in quotes) ")
@click.option(
    "--comment-char", 
    type=str, 
    default="#", 
    help="The first character of comment lines")
@click.option(
    "--send-header-to", 
    type=click.Choice(['dups', 'dedup', 'both', 'none']),
    default="both", 
    help="Which of the outputs should receive header and comment lines")
@click.option(
    "--c1", 
    type=int, 
    default=_pairsam_format.COL_C1,  
    help='Chrom 1 column; default {}'.format(_pairsam_format.COL_C1))
@click.option(
    "--c2", 
    type=int, 
    default=_pairsam_format.COL_C2,  
    help='Chrom 2 column; default {}'.format(_pairsam_format.COL_C2))
@click.option(
    "--p1", 
    type=int, 
    default=_pairsam_format.COL_P1,  
    help='Position 1 column; default {}'.format(_pairsam_format.COL_P1))
@click.option(
    "--p2", 
    type=int, 
    default=_pairsam_format.COL_P2,  
    help='Position 2 column; default {}'.format(_pairsam_format.COL_P2))
@click.option(
    "--s1", 
    type=int, 
    default=_pairsam_format.COL_S1,  
    help='Strand 1 column; default {}'.format(_pairsam_format.COL_S1))
@click.option(
    "--s2", 
    type=int, 
    default=_pairsam_format.COL_S2,  
    help='Strand 2 column; default {}'.format(_pairsam_format.COL_S2))
@click.option(
    "--unmapped-chrom", 
    type=str, 
    default=_pairsam_format.UNMAPPED_CHROM,  
    help='Placeholder for a chromosome on an unmapped side; default {}'.format(_pairsam_format.UNMAPPED_CHROM))

def dedup(pairsam_path, output, output_dups,
    stats_file,
    max_mismatch, method, 
    sep, comment_char, send_header_to,
    c1, c2, p1, p2, s1, s2, unmapped_chrom
    ):
    '''find and remove PCR duplicates.

    Find PCR duplicates in an upper-triangular flipped sorted pairs/pairsam 
    file. Allow for a +/-N bp mismatch at each side of duplicated molecules.

    PAIRSAM_PATH : input triu-flipped sorted .pairs or .pairsam file.  If the
    path ends with .gz, the input is gzip-decompressed. By default, the input 
    is read from stdin.
    '''
    dedup_py(pairsam_path, output, output_dups,
        stats_file,
        max_mismatch, method, 
        sep, comment_char, send_header_to,
        c1, c2, p1, p2, s1, s2, unmapped_chrom
        )


def dedup_py(pairsam_path, output, output_dups,
    stats_file,
    max_mismatch, method, 
    sep, comment_char, send_header_to,
    c1, c2, p1, p2, s1, s2, unmapped_chrom
    ):
    sep = ast.literal_eval('"""' + sep + '"""')
    send_header_to_dedup = send_header_to in ['both', 'dedup']
    send_header_to_dup = send_header_to in ['both', 'dups']

    instream = (_io.open_bgzip(pairsam_path, mode='r') 
                if pairsam_path else sys.stdin)
    outstream = (_io.open_bgzip(output, mode='w') 
                 if output else sys.stdout)
    outstream_dups = (_io.open_bgzip(output_dups, mode='w') 
                      if output_dups else None)

    header, body_stream = _headerops.get_header(instream)
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)

    if send_header_to_dedup:
        outstream.writelines((l+'\n' for l in header))
    if send_header_to_dup and outstream_dups:
        outstream_dups.writelines((l+'\n' for l in header))

    n_unmapped, n_dups, n_nodups = streaming_dedup(
        method, max_mismatch, sep, 
        c1, c2, p1, p2, s1, s2, unmapped_chrom,
        body_stream, outstream, outstream_dups)

    if stats_file:
        stat_f = _io.open_bgzip(stats_file, mode='a') 
        stat_f.write('{}\t{}\n'.format('n_unmapped', n_unmapped))
        stat_f.write('{}\t{}\n'.format('n_dups', n_dups))
        stat_f.write('{}\t{}\n'.format('n_nodups', n_nodups))
        stat_f.close()

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()

    if outstream_dups:
        outstream_dups.close()


def fetchadd(key, mydict):
    key = key.strip()
    if key not in mydict:
        mydict[key] = len(mydict)
    return mydict[key]


def ar(mylist, val):
    return np.array(mylist, dtype={8: np.int8, 16: np.int16, 32: np.int32}[val])
    

def streaming_dedup(
        method, max_mismatch, sep,
        c1ind, c2ind, p1ind, p2ind, s1ind, s2ind,
        unmapped_chrom,
        instream, outstream, outstream_dups):

    maxind = max(c1ind, c2ind, p1ind, p2ind, s1ind, s2ind)

    dd = _dedup.OnlineDuplicateDetector(method, max_mismatch, returnData=False)

    c1 = []; c2 = []; p1 = []; p2 = []; s1 = []; s2 = []
    lines = []
    chromDict = {}
    strandDict = {}
    n_unmapped = 0
    n_dups = 0
    n_nodups = 0

    while True: 
        line = next(instream, None)
        stripline = line.strip() if line else None

        if line:
            if not stripline: 
                warnings.warn("Empty line detected not at the end of the file")
                continue    

            words = line.split(sep)
            if len(words) <= maxind:
                raise ValueError(
                    "Error parsing line {}: ".format(line)
                    + " expected {} words, got {}".format(maxind, len(words)))
                
            if ((words[c1ind] == unmapped_chrom)
                or (words[c2ind] == unmapped_chrom)):

                outstream.write(line)  
                n_unmapped += 1
                    
            else:
                lines.append(line)
                c1.append(fetchadd(words[c1ind], chromDict))
                c2.append(fetchadd(words[c2ind], chromDict))
                p1.append(int(words[p1ind]))
                p2.append(int(words[p2ind]))
                s1.append(fetchadd(words[s1ind], strandDict))
                s2.append(fetchadd(words[s2ind], strandDict))

        if (not line) or (len(c1) == MAX_LEN):
            res = dd.push(ar(c1, 8), 
                          ar(c2, 8), 
                          ar(p1, 32), 
                          ar(p2, 32), 
                          ar(s1, 8), 
                          ar(s2, 8))
            if not line:
                res = np.concatenate([res, dd.finish()])

            for newline, remove in zip(lines[:len(res)], res):
                if not remove:
                    outstream.write(newline)  
                    n_nodups += 1
                else:
                    n_dups += 1
                    if outstream_dups:
                        outstream_dups.write(newline)
                    
            c1 = []; c2 = []; p1 = []; p2 = []; s1 = []; s2 = []
            lines = lines[len(res):]
            if not line:
                if(len(lines) != 0):                
                    raise ValueError(
                        "{} lines left in the buffer, ".format(len(lines))
                        + "should be none;"
                        + "something went terribly wrong")
                break

    return n_unmapped, n_dups, n_nodups


if __name__ == '__main__':
    dedup()
