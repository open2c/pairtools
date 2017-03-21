#!/usr/bin/env python
# -*- coding: utf-8  -*-
import argparse
import warnings
import sys
import ast 

import numpy as np
import pyximport; pyximport.install()
from _dedup import OnlineDuplicateDetector
from _distiller_common import open_bgzip, DISTILLER_VERSION, \
    append_pg_to_sam_header, get_header


# you don't need to load more than 10k lines at a time b/c you get out of the 
# CPU cache, so this parameter is not adjustable
MAX_LEN = 10000 


def main():
    parser = argparse.ArgumentParser(
        'Remove PCR duplicates from an upper-triangular flipped sorted '
        'pairs/pairsam file. Allow for a +/-N bp mismatch at each side of '
        'duplicated molecules.' )
    parser.add_argument(
        "--max-mismatch",
        type=int, 
        default=3,
        help='Pairs with both sides mapped within this distance (bp) from each '
             'other are considered duplicates.')
    parser.add_argument(
        '--sum',
        dest='METHOD', 
        action='store_const',
        const="sum", 
        default="max",  
        help='define the mismatch as the sum of the mismatches of the genomic '
            'locations of the both sides of the two compared molecules '
            '(default: use max of the two mismatches)')
    parser.add_argument(
        '--input',
        type=str, 
        default="",
        help='input triu-flipped sorted pairs or pairsam file.'
            ' If the path ends with .gz, the input is gzip-decompressed.'
            ' By default, the input is read from stdin.')
    parser.add_argument(
        "--output", 
        type=str, 
        default="", 
        help='output file for pairs after duplicate removal.'
            ' If the path ends with .gz, the output is bgzip-compressed.'
            ' By default, the output is printed into stdout.')
    parser.add_argument(
        "--output-dups",
        type=str, 
        default="", 
        help='output file for duplicates. '
            ' If the path ends with .gz, the output is bgzip-compressed.'
            ' By default, duplicates are dropped.'
            )
    parser.add_argument(
        "--sep",
        type=str, 
        default=r"\v", 
        help=r"Separator (\t, \v, etc. characters are "
              "supported, pass them in quotes) ")
    parser.add_argument(
        "--comment-char", 
        type=str, 
        default="#", 
        help="The first character of comment lines")
    parser.add_argument(
        "--send-header-to", 
        type=str, 
        default="both", 
        choices=['dups', 'dedup', 'both', 'none'], 
        help="Which of the outputs should receive header and comment lines")
    parser.add_argument(
        "--c1", 
        type=int, 
        default=0,  
        help='Chrom 1 column; default 0')
    parser.add_argument(
        "--c2", 
        type=int, 
        default=3,  
        help='Chrom 2 column; default 3')
    parser.add_argument(
        "--p1", 
        type=int, 
        default=1,  
        help='Position 1 column; default 1')
    parser.add_argument(
        "--p2", 
        type=int, 
        default=4,  
        help='Position 2 column; default 4')
    parser.add_argument(
        "--s1", 
        type=int, 
        default=2,  
        help='Strand 1 column; default 2')
    parser.add_argument(
        "--s2", 
        type=int, 
        default=5,  
        help='Strand 2 column; default 5')
    args = vars(parser.parse_args())

    sep = ast.literal_eval('"""' + args['sep'] + '"""')
    method = args['METHOD']
    max_mismatch = args['max_mismatch']
    comment_char = args['comment_char']
    send_header_to_dedup = args['send_header_to'] in ['both', 'dedup']
    send_header_to_dup = args['send_header_to'] in ['both', 'dups']
    c1ind = args['c1']
    c2ind = args['c2']
    p1ind = args['p1']
    p2ind = args['p2']
    s1ind = args['s1']
    s2ind = args['s2']

    instream = (open_bgzip(args['input'], mode='r') 
                if args['input'] else sys.stdin)
    outstream = (open_bgzip(args['output'], mode='w') 
                 if args['output'] else sys.stdout)
    outstream_dups = (open_bgzip(args['output_dups'], mode='w') 
                      if args['output_dups'] else None)

    header, pairsam_body_stream = get_header(instream)
    header = append_pg_to_sam_header(
        header,
        {'ID': 'pairs_dedup',
         'PN': 'pairs_dedup',
         'VN': DISTILLER_VERSION,
         'CL': ' '.join(sys.argv)
         })

    if send_header_to_dedup:
        outstream.writelines(header)
    if send_header_to_dup:
        outstream_dups.writelines(header)

    streaming_dedup(
        method, max_mismatch, sep, 
        c1ind, c2ind, p1ind, p2ind, s1ind, s2ind,
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
    main()
