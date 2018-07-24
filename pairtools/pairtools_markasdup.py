#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import pipes
import click

from . import _fileio, _pairsam_format, cli, _headerops, common_io_options

UTIL_NAME = 'pairtools_markasdup'

@cli.command()
@click.argument(
    'pairsam_path', 
    type=str,
    required=False)
@click.option(
    "-o", "--output", 
    type=str, 
    default="", 
    help='output .pairsam file.'
        ' If the path ends with .gz or .lz4, the output is pbgzip-/lz4c-compressed.'
        ' By default, the output is printed into stdout.')

@common_io_options

def markasdup(pairsam_path, output, **kwargs):
    '''Tag pairs as duplicates.

    Change the type of all pairs inside a .pairs/.pairsam file to DD. If sam
    entries are present, change the pair type in the Yt SAM tag to 'Yt:Z:DD'.

    PAIRSAM_PATH : input .pairs/.pairsam file. If the path ends with .gz, the 
    input is gzip-decompressed. By default, the input is read from stdin.
    '''
    markasdup_py(pairsam_path, output, **kwargs)

def markasdup_py(pairsam_path, output, **kwargs):
    instream = (_fileio.auto_open(pairsam_path, mode='r', 
                                  nproc=kwargs.get('nproc_in'),
                                  command=kwargs.get('cmd_in', None)) 
                if pairsam_path else sys.stdin)
    outstream = (_fileio.auto_open(output, mode='w', 
                                   nproc=kwargs.get('nproc_out'),
                                   command=kwargs.get('cmd_out', None)) 
                 if output else sys.stdout)

    header, body_stream = _headerops.get_header(instream)
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    outstream.writelines((l+'\n' for l in header))

    for line in body_stream:
        cols = line.rstrip().split(_pairsam_format.PAIRSAM_SEP)
        mark_split_pair_as_dup(cols)

        outstream.write(_pairsam_format.PAIRSAM_SEP.join(cols))
        outstream.write('\n')

    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()

def mark_split_pair_as_dup(cols):
    # if the original columns ended with a new line, the marked columns
    # should as well.
    original_has_newline = cols[-1].endswith('\n')

    cols[_pairsam_format.COL_PTYPE] = 'DD'
    
    if (len(cols) > _pairsam_format.COL_SAM1) and (len(cols) > _pairsam_format.COL_SAM2):
        for i in (_pairsam_format.COL_SAM1,
                  _pairsam_format.COL_SAM2):
                
            # split each sam column into sam entries, tag and assemble back
            cols[i] = _pairsam_format.INTER_SAM_SEP.join(
                [mark_sam_as_dup(sam) 
                 for sam in cols[i].split(_pairsam_format.INTER_SAM_SEP)
                ])
                
    if original_has_newline and not cols[-1].endswith('\n'):
        cols[-1] = cols[-1]+'\n'
    return cols

def mark_sam_as_dup(sam):
    '''Tag the binary flag and the optional pair type field of a sam entry
    as a PCR duplicate.'''
    samcols = sam.split(_pairsam_format.SAM_SEP)

    if len(samcols) == 1:
        return sam

    samcols[1] = str(int(samcols[1]) | 1024)

    for j in range(11, len(samcols)):
        if samcols[j].startswith('Yt:Z:'):
            samcols[j] = 'Yt:Z:DD'
    return _pairsam_format.SAM_SEP.join(samcols)


if __name__ == '__main__':
    markasdup()
