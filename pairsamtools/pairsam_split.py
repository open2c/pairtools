#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import pipes
import click

from . import _fileio, _pairsam_format, _headerops, cli, common_io_options

UTIL_NAME = 'pairsam_split'

@cli.command()
@click.argument(
    'pairsam_path', 
    type=str,
    required=False)

@click.option(
    "--output-pairs", 
    type=str, 
    default="", 
    help='output pairs file.'
        ' If the path ends with .gz or .lz4, the output is pbgzip-/lz4c-compressed.'
        ' If -, pairs are printed to stdout.'
        ' If not specified, pairs are dropped.')
@click.option(
    "--output-sam", 
    type=str, 
    default="", 
    help='output sam file.'
        ' If the path ends with .bam, the output is compressed into a bam file.'
        ' If -, sam entries are printed to stdout.'
        ' If not specified, sam entries are dropped.')

@common_io_options

def split(pairsam_path, output_pairs, output_sam, **kwargs):
    '''split a .pairsam file into pairs and sam.

    PAIRSAM_PATH : input .pairsam file. If the path ends with .gz or .lz4, the
    input is decompressed by pbgzip or lz4c. By default, the input is read from 
    stdin.
    '''
    split_py(pairsam_path, output_pairs, output_sam, **kwargs)


def split_py(pairsam_path, output_pairs, output_sam, **kwargs):
    instream = (_fileio.auto_open(pairsam_path, mode='r', 
                                  nproc=kwargs.get('nproc_in'),
                                  command=kwargs.get('cmd_in', None)) 
                if pairsam_path else sys.stdin)

    # Output streams
    if (not output_pairs) and (not output_sam):
        raise ValueError('At least one output (pairs and/or sam) must be specified!')
    if (output_pairs == '-') and (output_sam == '-'):
        raise ValueError('Only one output (pairs or sam) can be printed in stdout!')

    outstream_pairs = (sys.stdout if (output_pairs=='-')
                       else (_fileio.auto_open(output_pairs, mode='w', 
                                               nproc=kwargs.get('nproc_out'),
                                               command=kwargs.get('cmd_out', None)) 
                             if output_pairs else None))
    outstream_sam = (sys.stdout if (output_sam=='-')
                     else (_fileio.auto_open(output_sam, mode='w',
                                             nproc=kwargs.get('nproc_out'),
                                             command=kwargs.get('cmd_out', None))
                           if output_sam else None))

    header, body_stream = _headerops.get_header(instream)
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)

    if outstream_pairs:
        outstream_pairs.writelines((l+'\n' for l in header))
    if outstream_sam:
        outstream_sam.writelines(
            (l[11:].strip()+'\n' for l in header if l.startswith('#samheader:')))

    # Split
    for line in body_stream:
        cols = line.rstrip().split(_pairsam_format.PAIRSAM_SEP)
        if outstream_pairs:
            # hard-coded tab separator to follow the DCIC pairs standard
            outstream_pairs.write('\t'.join(cols[:_pairsam_format.COL_SAM1]))
            outstream_pairs.write('\n')
        
        if (outstream_sam 
            and (len(cols) > _pairsam_format.COL_SAM1) 
            and (len(cols) > _pairsam_format.COL_SAM2)):

            for col in (cols[_pairsam_format.COL_SAM1],
                        cols[_pairsam_format.COL_SAM2]):
                if col != '.':
                    for sam_entry in col.split(_pairsam_format.INTER_SAM_SEP):
                        outstream_sam.write(sam_entry.replace(_pairsam_format.SAM_SEP,'\t'))
                        outstream_sam.write('\n')

    if outstream_pairs and outstream_pairs != sys.stdout:
        outstream_pairs.close()

    if outstream_sam and outstream_sam != sys.stdout:
        outstream_sam.close()


if __name__ == '__main__':
    split()
