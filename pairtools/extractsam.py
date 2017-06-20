#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import pipes
import click

from . import _io, _pairsam_format, _headerops, cli

UTIL_NAME = 'pairtools_extractsam'

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
        ' If the path ends with .gz, the output is bgzip-compressed.'
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

@click.option(
    "--nproc", 
    type=int, 
    default=8, 
    show_default=True,
    help='Number of processes to split the work between.'
    )

def extractsam(pairsam_path, output_pairs, output_sam, nproc):
    '''Extract SAM entries from a .pairsam file a form a .pairs file and a 
    .sam file.

    PAIRSAM_PATH : input .pairsam file. If the path ends with .gz, the input is
    gzip-decompressed. By default, the input is read from stdin.
    '''
    extractsam_py(pairsam_path, output_pairs, output_sam, nproc)


def extractsam_py(pairsam_path, output_pairs, output_sam, nproc):
    instream = (_io.open_bgzip(pairsam_path, mode='r', nproc=nproc) 
                if pairsam_path else sys.stdin)

    # Output streams
    if (not output_pairs) and (not output_sam):
        raise Exception('At least one output (pairs and/or sam) must be specified!')
    if (output_pairs == '-') and (output_sam == '-'):
        raise Exception('Only one output (pairs or sam) can be printed in stdout!')

    outstream_pairs = (sys.stdout if (output_pairs=='-')
                  else _io.open_bgzip(output_pairs, mode='w', nproc=nproc) if output_pairs
                  else None)
    outstream_sam = (sys.stdout if (output_sam=='-')
                else _io.open_sam_or_bam(output_sam, mode='w', nproc=nproc) if output_sam
                else None)

    header, body_stream = _headerops.get_header(instream)
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)

    if outstream_pairs:
        outstream_pairs.writelines((l+'\n' for l in header))
    if outstream_sam:
        outstream_sam.writelines(
            (l[11:].strip()+'\n' for l in header if l.startswith('#samheader:')))

    # Split
    for line in body_stream:
        cols = line[:-1].split(_pairsam_format.PAIRSAM_SEP)
        if outstream_pairs:
            # hard-coded tab separator to follow the DCIC pairs standard
            outstream_pairs.write('\t'.join(cols[:_pairsam_format.COL_SAM1] 
                                            + cols[_pairsam_format.COL_SAM2+1:]))
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
    extractsam()
