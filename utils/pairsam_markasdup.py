#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import pipes
import sys

import _distiller_common

UTIL_NAME = 'pairsam_markasdup'


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

    args = vars(parser.parse_args())
    
    instream = (_distiller_common.open_bgzip(args['input'], mode='r') 
                if args['input'] else sys.stdin)
    outstream = (_distiller_common.open_bgzip(args['output'], mode='w') 
                 if args['output'] else sys.stdout)
 
    header, pairsam_body_stream = _distiller_common.get_header(instream)
    header = _distiller_common.append_pg_to_sam_header(
        header,
        {'ID': UTIL_NAME,
         'PN': UTIL_NAME,
         'VN': _distiller_common.DISTILLER_VERSION,
         'CL': ' '.join(sys.argv)
         })

    outstream.writelines(header)

    for line in pairsam_body_stream:
        cols = line[:-1].split('\v')
        cols[_distiller_common.COL_PTYPE] = 'DD'
        
        for i in (_distiller_common.COL_SAM1,
                  _distiller_common.COL_SAM2):
                
            # split each sam column into sam entries, tag and assemble back
            cols[i] = _distiller_common.SAM_ENTRY_SEP.join(
                [mark_sam_as_dup(sam) 
                 for sam in cols[i].split(_distiller_common.SAM_ENTRY_SEP)
                ])

        outstream.write('\v'.join(cols))
        outstream.write('\n')

    if hasattr(instream, 'close'):
        instream.close()
    if hasattr(outstream, 'close'):
        outstream.close()

def mark_sam_as_dup(sam):
    '''Tag the binary flag and the optional pair type field of a sam entry
    as a PCR duplicate.'''
    samcols = sam.split('\t')
    samcols[1] = str(int(samcols[1]) | 1024)

    for j in range(11, len(samcols)):
        if samcols[j].startswith('Yt:Z:'):
            samcols[j] = 'Yt:Z:DD'
    return '\t'.join(samcols)


if __name__ == '__main__':
    main()
