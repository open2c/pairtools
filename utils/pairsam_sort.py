#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import sys
import argparse
import subprocess

from _distiller_common import open_bgzip, DISTILLER_VERSION, \
    append_pg_to_sam_header, get_header

def main():
    parser = argparse.ArgumentParser(
        description='''Sort a pairsam file. The resulting order is lexicographic
        along chrom1 and chrom2, numeric along pos1 and pos2 and lexicographic
        along pair_type.
        '''
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
    
    instream = (open_bgzip(args['input'], mode='r') 
                if args['input'] else sys.stdin)
    outstream = (open_bgzip(args['output'], mode='w') 
                 if args['output'] else sys.stdout)

    header, pairsam_body_stream = get_header(instream)
    header = append_pg_to_sam_header(
        header,
        {'ID': 'pairsam_sort',
         'PN': 'pairsam_sort',
         'VN': DISTILLER_VERSION,
         'CL': ' '.join(sys.argv)
         })

    outstream.writelines(header)
 
    if hasattr(outstream, 'close'):
        outstream.close()

    command = r'''
        /bin/bash -c 'sort -k 1,1 -k 4,4 -k 2,2n -k 5,5n -k 8,8 
        --field-separator=$'\''\v'\'' 
        '''.replace('\n','')
    if args['output'].endswith('.gz'):
        command += '| bgzip -c'
    if args['output']:
        command += ' >> ' + args['output']
    command += "'"

    with subprocess.Popen(command, stdin=subprocess.PIPE, bufsize=-1, shell=True) as process:
        stdin_wrapper = io.TextIOWrapper(process.stdin, 'utf-8')
        for line in pairsam_body_stream:
            stdin_wrapper.write(line)
        stdin_wrapper.flush()

    if hasattr(instream, 'close'):
        instream.close()


if __name__ == '__main__':
    main()
