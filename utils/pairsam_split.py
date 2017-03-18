#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import pipes
import sys
from _distiller_common import open_bam_or_sam, open_bgzip


def main():
    parser = argparse.ArgumentParser(
        'Splits a .pairsam file into pairs and sam entries')
    parser.add_argument(
        'infile', nargs='?', 
        type=argparse.FileType('r'), 
        default=sys.stdin)
    parser.add_argument(
        "--out-pairs", 
        type=str, 
        required=True)
    parser.add_argument(
        "--out-sam", 
        type=str, 
        required=True)
    parser.add_argument(
        "--comment-char", 
        type=str, 
        default="#", 
        help="The first character of comment lines")
    args = vars(parser.parse_args())

    # Output streams
    if args['out_pairs'] is None:
        pairs_file = sys.stdout
    else:
        pairs_file = open_bgzip(args['out_pairs'], mode='w')
    sam_file = open_sam_or_bam(args['out_sam'], 'w')

    # Input pairsam
    instream = args['infile']
    comment_char = args['comment_char']

    # Split
    for line in instream.readlines():
        if line.startswith(comment_char):
            if line.startswith(comment_char+'@'):
                sam_file.write(line[len(comment_char):])
            else:
                pairs_file.write(line)
            continue

        cols = line[:-1].split('\v')
        pairs_file.write('\t'.join(cols[:8]))
        pairs_file.write('\n')
        
        for col in cols[8:]:
            sam_file.write(col)
            sam_file.write('\n')

    if hasattr(pairs_file, 'close'):
        pairs_file.close()

    if hasattr(sam_file, 'close'):
        sam_file.close()


if __name__ == '__main__':
    main()
