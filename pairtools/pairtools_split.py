#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import pipes
import click

from . import _fileio, _pairsam_format, _headerops, cli, common_io_options

UTIL_NAME = 'pairtools_split'

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
    '''Split a .pairsam file into .pairs and .sam.

    Restore a .sam file from sam1 and sam2 fields of a .pairsam file. Create 
    a .pairs file without sam1/sam2 fields.

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
    columns = _headerops.extract_column_names(header)

    has_sams = False
    if columns:
        # trust the column order specified in the header
        if ('sam1' in columns) and ('sam2' in columns):
            sam1col = columns.index('sam1') 
            sam2col = columns.index('sam2')
            columns.pop(max(sam1col, sam2col))
            columns.pop(min(sam1col, sam2col))
            header = _headerops._update_header_entry(
                header, 'columns', ' '.join(columns))
            has_sams = True
        elif ('sam1' in columns) != ('sam1' in columns):
            raise ValueError(
                'According to the #columns header field only one sam entry is present')
    else:
        # assume that the file has sam columns and follows the pairsam format
        sam1col = _pairsam_format.COL_SAM1
        sam2col = _pairsam_format.COL_SAM2
        has_sams = True

    if outstream_pairs:
        outstream_pairs.writelines((l+'\n' for l in header))
    if outstream_sam:
        outstream_sam.writelines(
            (l[11:].strip()+'\n' for l in header if l.startswith('#samheader:')))

    # Split
    sam1 = None
    sam2 = None
    for line in body_stream:
        cols = line.rstrip().split(_pairsam_format.PAIRSAM_SEP)
        if has_sams:
            if sam1col < sam2col:
                sam2 = cols.pop(sam2col)
                sam1 = cols.pop(sam1col)
            else:
                sam1 = cols.pop(sam1col)
                sam2 = cols.pop(sam2col)

        if outstream_pairs:
            # hard-coded tab separator to follow the DCIC pairs standard
            outstream_pairs.write('\t'.join(cols))
            outstream_pairs.write('\n')
        
        if (outstream_sam and has_sams):
            for col in (sam1, sam2):
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
