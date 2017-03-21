#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import pipes
import sys

from _distiller_common import open_bgzip, DISTILLER_VERSION, \
    append_pg_to_sam_header, get_header

def main():
    parser = argparse.ArgumentParser(
        description='Tags every line of a pairsam with a duplicate tag'
    )
    parser.add_argument(
        '--input',
        type=str, 
        default="",
        help='input pairsam file.'
            ' If the path ends with .gz, the input is gzip-decompressed.'
            ' By default, the input is read from stdin.')

    parser.add_argument(
        "--output", 
        type=str, 
        default="", 
        help='output pairsam file.'
            ' If the path ends with .gz, the output is bgzip-compressed.'
            ' By default, the output is printed into stdout.')

    parser.add_argument(
        "--comment-char", 
        type=str, 
        default="#",
        help="The first character of comment lines")
    args = vars(parser.parse_args())
    
    instream = (open_bgzip(args['input'], mode='r') 
                if args['input'] else sys.stdin)
    outstream = (open_bgzip(args['output'], mode='w') 
                 if args['output'] else sys.stdout)
 
    comment_char = args['comment_char']

    header, pairsam_body_stream = get_header(instream)
    header = append_pg_to_sam_header(
        header,
        {'ID': 'pairsam_markasdup',
         'PN': 'pairsam_markasdup',
         'VN': DISTILLER_VERSION,
         'CL': ' '.join(sys.argv)
         })

    outstream.writelines(header)

    for line in pairsam_body_stream:
        cols = line[:-1].split('\v')
        cols[7] = 'DD'
        
        for i in range(8, len(cols)):
            sam = cols[i]
            samcols = sam.split('\t')
            samcols[1] = str(int(samcols[1]) | 1024)

            for j in range(11, len(samcols)):
                if samcols[j].startswith('Yt:Z:'):
                    samcols[j] = 'Yt:Z:DD'
            
            cols[i] = '\t'.join(samcols)

        outstream.write('\v'.join(cols))
        outstream.write('\n')


if __name__ == '__main__':
    main()
