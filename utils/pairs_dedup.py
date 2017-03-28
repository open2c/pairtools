#!/usr/bin/env python
# -*- coding: utf-8  -*-
import sys
import ast 
import warnings

import click

import numpy as np
import pyximport; pyximport.install()
from _dedup import OnlineDuplicateDetector

import _distiller_common

UTIL_NAME = 'pairs_dedup'


# you don't need to load more than 10k lines at a time b/c you get out of the 
# CPU cache, so this parameter is not adjustable
MAX_LEN = 10000 


@click.command()
@click.option(
    '--input',
    type=str, 
    default="",
    help='input triu-flipped sorted pairs or pairsam file.'
        ' If the path ends with .gz, the input is gzip-decompressed.'
        ' By default, the input is read from stdin.')
@click.option(
    "--output", 
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
        ' By default, duplicates are dropped.'
        )
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
    show_default=True,
        )
@click.option(
    "--sep",
    type=str, 
    default=r"\v", 
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
    default=_distiller_common.COL_C1,  
    help='Chrom 1 column; default {}'.format(_distiller_common.COL_C1))
@click.option(
    "--c2", 
    type=int, 
    default=_distiller_common.COL_C2,  
    help='Chrom 2 column; default {}'.format(_distiller_common.COL_C2))
@click.option(
    "--p1", 
    type=int, 
    default=_distiller_common.COL_P1,  
    help='Position 1 column; default {}'.format(_distiller_common.COL_P1))
@click.option(
    "--p2", 
    type=int, 
    default=_distiller_common.COL_P2,  
    help='Position 2 column; default {}'.format(_distiller_common.COL_P2))
@click.option(
    "--s1", 
    type=int, 
    default=_distiller_common.COL_S1,  
    help='Strand 1 column; default {}'.format(_distiller_common.COL_S1))
@click.option(
    "--s2", 
    type=int, 
    default=_distiller_common.COL_S2,  
    help='Strand 2 column; default {}'.format(_distiller_common.COL_S2))

def dedup(input, output, output_dups, max_mismatch, method, 
    sep, comment_char, send_header_to,
    c1, c2, p1, p2, s1, s2
    ):
    '''Remove PCR duplicates from an upper-triangular flipped sorted 
    pairs/pairsam file. Allow for a +/-N bp mismatch at each side of 
    duplicated molecules.'''

    sep = ast.literal_eval('"""' + sep + '"""')
    send_header_to_dedup = send_header_to in ['both', 'dedup']
    send_header_to_dup = send_header_to in ['both', 'dups']

    instream = (_distiller_common.open_bgzip(input, mode='r') 
                if input else sys.stdin)
    outstream = (_distiller_common.open_bgzip(output, mode='w') 
                 if output else sys.stdout)
    outstream_dups = (_distiller_common.open_bgzip(output_dups, mode='w') 
                      if output_dups else None)

    header, pairsam_body_stream = _distiller_common.get_header(instream)
    header = _distiller_common.append_pg_to_sam_header(
        header,
        {'ID': UTIL_NAME,
         'PN': UTIL_NAME,
         'VN': _distiller_common.DISTILLER_VERSION,
         'CL': ' '.join(sys.argv)
         })

    if send_header_to_dedup:
        outstream.writelines(header)
    if send_header_to_dup and outstream_dups:
        outstream_dups.writelines(header)

    streaming_dedup(
        method, max_mismatch, sep, 
        c1, c2, p1, p2, s1, s2,
        pairsam_body_stream, outstream, outstream_dups)

    if hasattr(instream, 'close'):
        instream.close()
    if hasattr(outstream, 'close'):
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
        instream, outstream, outstream_dups):

    maxind = max(c1ind, c2ind, p1ind, p2ind, s1ind, s2ind)

    dd = OnlineDuplicateDetector(method, max_mismatch, returnData=False)

    c1 = []; c2 = []; p1 = []; p2 = []; s1 = []; s2 = []
    lines = []
    chromDict = {}
    strandDict = {}

    while True: 
        line = next(instream, None)
        stripline = line.strip() if line else None

        if line:
            if not stripline: 
                warnings.warn("Empty line detected not at the end of the file")
                continue    

            lines.append(line)
            words = line.split(sep)
            if len(words) <= maxind:
                raise ValueError(
                    "Error parsing line {}: ".format(line)
                    + " expected {} words, got {}".format(maxind, len(words)))
                
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
                else:
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


if __name__ == '__main__':
    dedup()
