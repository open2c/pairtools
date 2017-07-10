#!/usr/bin/env python
import sys
import glob
import math
import subprocess
import click

from . import _io, _pairsam_format, _headerops, cli

UTIL_NAME = 'pairsam_merge'


@cli.command()
@click.argument(
    'pairsam_path', 
    nargs=-1, 
    type=str,
    )
@click.option(
    "-o", "--output", 
    type=str, 
    default="", 
    help='output file.'
        ' If the path ends with .gz, the output is bgzip-compressed.'
        ' By default, the output is printed into stdout.')
@click.option(
    "--npzipin", 
    type=int, 
    default=2, 
    help='Number of threads per pbgzip instance in the input stream (decompression).'
        ' Decompression of all individual chunks starts simultaneously,'
        ' so, take into account total number of merging chunks to adjust this flag',
    show_default=True,
    )
@click.option(
    "--npzipout", 
    type=int, 
    default=4, 
    help='Number of threads per pbgzip instance in the output stream (compression).',
    show_default=True,
    )
@click.option(
    "--npmerge", 
    type=int, 
    default=2, 
    help='Number of threads per sort --merge.',
    show_default=True,
    )
@click.option(
    "--batchmerge", 
    type=int, 
    default=0, 
    help='Merging batch size for sort --merge --batch-size'
        ' Default value implies all chunks are merged together simultaneously in a single batch (no batching)'
        ' This flag controls the batch size only.'
        ' It does not control the number of batches being merged simultaneously'
        ' This flag does not control the number of input decompression streams in any way.',
    show_default=True,
    )
@click.option(
    "--tempdir", 
    type=str, 
    default=".", 
    help='Path to store temporary files that sort --merge generates'
        ' Use it only with batchmerge flag.',
    )



def merge(pairsam_path, output, npzipin, npzipout, npmerge, batchmerge, tempdir):
    """merge sorted pairs/pairsam files. 

    Merge triu-flipped sorted pairs/pairsam files. If present, the @SQ records 
    of the SAM header must be identical; the sorting order of 
    these lines is taken from the first file in the list. 
    The ID fields of the @PG records of the SAM header are modified with a
    numeric suffix to produce unique records.
    The other unique SAM and non-SAM header lines are copied into the output header.

    PAIRSAM_PATH : upper-triangular flipped sorted pairs/pairsam files to merge
    or a group/groups of .pairsam files specified by a wildcard
    
    """
    merge_py(pairsam_path, output, npzipin, npzipout, npmerge, batchmerge, tempdir)


def merge_py(pairsam_path, output, npzipin, npzipout, npmerge, batchmerge, tempdir):
    paths = sum([glob.glob(mask) for mask in pairsam_path], [])

    outstream = (_io.open_bgzip(output, mode='w', nproc=npzipout) 
                 if output else sys.stdout)

    headers = []
    for path in paths:
        f = _io.open_bgzip(path, mode='r')
        h, _ = _headerops.get_header(f)
        headers.append(h)
        f.close()
    merged_header = _headerops.merge_headers(headers)
    merged_header = _headerops.append_new_pg(
        merged_header, ID=UTIL_NAME, PN=UTIL_NAME)

    outstream.writelines((l+'\n' for l in merged_header))
    outstream.flush()
 
    command = r'''
        /bin/bash -c 'sort -k {0},{0} -k {1},{1} -k {2},{2}n -k {3},{3}n -k {4},{4} 
        --merge  --field-separator=$'\''{5}'\''
        {6}
        {7}
        {8}
        -S 1G
        '''.replace('\n',' ').format(
                _pairsam_format.COL_C1+1, 
                _pairsam_format.COL_C2+1, 
                _pairsam_format.COL_P1+1, 
                _pairsam_format.COL_P2+1,
                _pairsam_format.COL_PTYPE+1,
                _pairsam_format.PAIRSAM_SEP_ESCAPE,
                ' --parallel={} '.format(npmerge) if npmerge > 1 else ' ',
                ' --temporary-directory={}'.format(tempdir),
                ' --batch-size={}'.format(batchmerge) if batchmerge > 1 else ' ',
                )
    for path in paths:
        if path.endswith('.gz'):
            command += r''' <(pbgzip -dc -n {} {} | sed -n -e '\''/^[^#]/,$p'\'')'''.format(npzipin, path)
        else:
            command += r''' <(sed -n -e '\''/^[^#]/,$p'\'' {})'''.format(path)
    command += "'"
    subprocess.call(command, shell=True, stdout=outstream)

    if outstream != sys.stdout:
        outstream.close()


if __name__ == '__main__':
    merge()
