#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import sys
import click
import subprocess

from . import _common, cli

UTIL_NAME = 'pairsam_sort'

@cli.command()

@click.option(
    '--input',
    type=str, 
    default="",
    help='input pairsam file.'
        ' If the path ends with .gz, the input is gzip-decompressed.'
        ' By default, the input is read from stdin.')

@click.option(
    "--output", 
    type=str, 
    default="", 
    help='output pairsam file.'
        ' If the path ends with .gz, the output is bgzip-compressed.'
        ' By default, the output is printed into stdout.')

def sort(input, output):
    '''Sort a pairsam file. The resulting order is lexicographic
    along chrom1 and chrom2, numeric along pos1 and pos2 and lexicographic
    along pair_type.
    '''

    instream = (_common.open_bgzip(input, mode='r') 
                if input else sys.stdin)
    outstream = (_common.open_bgzip(output, mode='w') 
                 if output else sys.stdout)

    header, pairsam_body_stream = _common.get_header(instream)
    header = _common.append_pg_to_sam_header(
        header,
        {'ID': UTIL_NAME,
         'PN': UTIL_NAME,
         'VN': _common.DISTILLER_VERSION,
         'CL': ' '.join(sys.argv)
         })

    outstream.writelines(header)
    outstream.flush()

    command = r'''
        /bin/bash -c 'sort 
        -k {0},{0} -k {1},{1} -k {2},{2}n -k {3},{3}n -k {4},{4} 
        --field-separator=$'\''\v'\'' 
        '''.replace('\n',' ').format(
                _common.COL_C1+1, 
                _common.COL_C2+1, 
                _common.COL_P1+1, 
                _common.COL_P2+1,
                _common.COL_PTYPE+1,
        )
    command += "'"

    with subprocess.Popen(
            command, stdin=subprocess.PIPE, bufsize=-1, shell=True,
            stdout=outstream) as process:
        stdin_wrapper = io.TextIOWrapper(process.stdin, 'utf-8')
        for line in pairsam_body_stream:
            stdin_wrapper.write(line)
        stdin_wrapper.flush()
        process.communicate()

    if hasattr(instream, 'close'):
        instream.close()


if __name__ == '__main__':
    sort()
