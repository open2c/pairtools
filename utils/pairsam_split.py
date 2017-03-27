#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import pipes
import sys

import _distiller_common

def main():
    parser = argparse.ArgumentParser(
        description='Splits a .pairsam file into pairs and sam entries')
    parser.add_argument(
        '--input',
        type=str, 
        default="",
        help='input pairsam file.'
            ' If the path ends with .gz, the input is gzip-decompressed.'
            ' By default, the input is read from stdin.')
    parser.add_argument(
        "--output-pairs", 
        type=str, 
        required=True)
    parser.add_argument(
        "--output-sam", 
        type=str, 
        required=True)
    parser.add_argument(
        "--comment-char", 
        type=str, 
        default="#", 
        help="The first character of comment lines")
    args = vars(parser.parse_args())

    instream = (_distiller_common.open_bgzip(args['input'], mode='r') 
                if args['input'] else sys.stdin)

    # Output streams
    pairs_file = _distiller_common.open_bgzip(args['output_pairs'], mode='w') 
    sam_file = _distiller_common.open_sam_or_bam(args['output_sam'], 'w')

    # Input pairsam
    comment_char = args['comment_char']

    # Split
    for line in instream.readlines():
        if line.startswith(comment_char):
            if line.startswith(comment_char+'@'):
                sam_file.write(line[len(comment_char):])

            pairs_file.write(line)
            continue

        cols = line[:-1].split('\v')
        pairs_file.write('\t'.join(cols[:_distiller_common.COL_SAM]))
        pairs_file.write('\n')
        
        for col in cols[_distiller_common.COL_SAM:]:
            sam_file.write(col)
            sam_file.write('\n')

    if hasattr(pairs_file, 'close'):
        pairs_file.close()

    if hasattr(sam_file, 'close'):
        sam_file.close()


if __name__ == '__main__':
    main()
