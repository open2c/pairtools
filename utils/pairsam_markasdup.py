#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import pipes
import sys


def main():
    parser = argparse.ArgumentParser(
        description='Tags every line of a pairsam with a duplicate tag'
    )
    parser.add_argument(
        'infile', nargs='?', 
        type=argparse.FileType('r'), 
        default=sys.stdin)
    parser.add_argument(
        "--comment-char", 
        type=str, 
        default="#",
        help="The first character of comment lines")
    args = vars(parser.parse_args())


    comment_char = args['comment_char']
    instream = args['infile']
    outstream = sys.stdout

    streaming_tag_dups(instream, outstream, comment_char)


def streaming_tag_dups(instream, outstream, comment_char):
    for line in instream.readlines():

        if line.startswith(comment_char):
            outstream.write(line)
            continue
        
        else:
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
