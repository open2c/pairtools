#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import pipes
import click

import _distiller_common

@click.command()
@click.option(
    '--input',
    type=str, 
    default="",
    help='input pairsam file.'
        ' If the path ends with .gz, the input is gzip-decompressed.'
        ' By default, the input is read from stdin.')
@click.argument(
    "output_pairs", 
    metavar='OUTPUT_PAIRS', 
    type=str, 
    )
@click.argument(
    "output_sam", 
    metavar='OUTPUT_SAM', 
    type=str, 
    )

def split(input, output_pairs, output_sam):
    '''Splits a .pairsam file into pairs and sam entries

    OUTPUT_PAIRS : output pairs file. If the path ends with .gz, the output is 
    bgzip-compressed.

    OUTPUT_SAM : output pairs file. If the path ends with .bam, the output is
    compressed into a bam file.

    '''
    instream = (_distiller_common.open_bgzip(input, mode='r') 
                if input else sys.stdin)

    # Output streams
    pairs_file = _distiller_common.open_bgzip(output_pairs, mode='w') 
    sam_file = _distiller_common.open_sam_or_bam(output_sam, 'w')


    # Split
    for line in instream.readlines():
        if line.startswith('#'):
            if line.startswith('#'+'@'):
                sam_file.write(line[len('#'):])

            pairs_file.write(line)
            continue

        cols = line[:-1].split('\v')
        pairs_file.write('\t'.join(cols[:_distiller_common.COL_SAM1]))
        pairs_file.write('\n')
        
        for col in (cols[_distiller_common.COL_SAM1],
                    cols[_distiller_common.COL_SAM2]):
            for sam_entry in col.split(_distiller_common.SAM_ENTRY_SEP):
                sam_file.write(sam_entry)
                sam_file.write('\n')

    if hasattr(pairs_file, 'close'):
        pairs_file.close()

    if hasattr(sam_file, 'close'):
        sam_file.close()


if __name__ == '__main__':
    split()
