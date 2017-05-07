#!/usr/bin/env python
import sys
import glob
import math
import subprocess
import click

from . import _io, _common, _headerops, cli

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
    "--nproc", 
    type=int, 
    default=8, 
    help='Number of threads per process.',
    show_default=True,
    )


def merge(pairsam_path, output, nproc):
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
    merge_py(pairsam_path, output, nproc)

def merge_py(pairsam_path, output, nproc):
    paths = sum([glob.glob(mask) for mask in pairsam_path], [])

    outstream = (_io.open_bgzip(output, mode='w', nproc=nproc) 
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
        -S 1G
        '''.replace('\n',' ').format(
                _common.COL_C1+1, 
                _common.COL_C2+1, 
                _common.COL_P1+1, 
                _common.COL_P2+1,
                _common.COL_PTYPE+1,
                _common.PAIRSAM_SEP_ESCAPE,
                ' --parallel={} '.format(nproc) if nproc > 1 else ' ',
                )
    for path in paths:
        if path.endswith('.gz'):
            command += r''' <(pbgzip -dc -n {} {} | sed -n -e '\''/^[^#]/,$p'\'')'''.format(nproc, path)
        else:
            command += r''' <(sed -n -e '\''/^[^#]/,$p'\'' {})'''.format(path)
    command += "'"
    subprocess.call(command, shell=True, stdout=outstream)

    if outstream != sys.stdout:
        outstream.close()


if __name__ == '__main__':
    merge()
