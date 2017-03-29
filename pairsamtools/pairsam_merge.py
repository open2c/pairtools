#!/usr/bin/env python
import sys
import glob
import subprocess
import click

from . import _common, cli

UTIL_NAME = 'pairsam_merge'

@cli.command()
@click.argument(
    'infile', 
    metavar='INFILE',
    nargs=-1, 
    type=str,
    )
@click.option(
    "--output", 
    type=str, 
    default="", 
    help='output file.'
        ' If the path ends with .gz, the output is bgzip-compressed.'
        ' By default, the output is printed into stdout.')

def merge(infile,output):
    """Merge multiple sorted pairsam files. 
    The @SQ records of the SAM header must be identical; the sorting order of 
    these lines is taken from the first file in the list. 
    The ID fields of the @PG records of the SAM header are modified with a
    numeric suffix to produce unique records.
    The other unique SAM and non-SAM header lines are copied into the output header.

    INFILE : a file to merge or a group of files specified by a wildcard
    
    """

    outstream = (_common.open_bgzip(output, mode='w') 
                 if output else sys.stdout)

    paths = sum([glob.glob(mask) for mask in infile], [])
    merged_header = form_merged_header(paths)

    merged_header = _common.append_pg_to_sam_header(
        merged_header,
        {'ID': UTIL_NAME,
         'PN': UTIL_NAME,
         'VN': _common.DISTILLER_VERSION,
         'CL': ' '.join(sys.argv)
         })

    outstream.writelines(merged_header)
    outstream.flush()
 
    command = r'''
        /bin/bash -c 'sort -k {0},{0} -k {1},{1} -k {2},{2}n -k {3},{3}n -k {4},{4} 
        --merge --field-separator=$'\''\v'\'' 
        '''.replace('\n',' ').format(
                _common.COL_C1+1, 
                _common.COL_C2+1, 
                _common.COL_P1+1, 
                _common.COL_P2+1,
                _common.COL_PTYPE+1,
                )
    for path in paths:
        if path.endswith('.gz'):
            command += r''' <(zcat {} | sed -n -e '\''/^[^#]/,$p'\'')'''.format(path)
        else:
            command += r''' <(sed -n -e '\''/^[^#]/,$p'\'' {})'''.format(path)
    command += "'"
    subprocess.call(command, shell=True, stdout=outstream)


def form_merged_header(paths):
    headers = []
    for path in paths:
        header = []
        f = _common.open_bgzip(path, mode='r')
        for line in f.readlines():
            if line and not line.isspace():
                if line.strip().startswith('#'):
                    header.append(line.strip())
                else:
                    break
        f.close()
        headers.append(header)

    # HD headers contain information that becomes invalid after processing
    # with distiller. Do not print into the output.
    HD_headers = [set(line for line in header if line.startswith('#@HD'))
                  for header in headers]

    SQ_headers = [set(line for line in header if line.startswith('#@SQ'))
                  for header in headers]
    common_sq_header = set.intersection(*SQ_headers)
    sq_headers_same = all([len(header) == len(common_sq_header) 
                           for header in SQ_headers])

    if not sq_headers_same:
        raise Exception('The SQ (sequence) lines of the sam headers are not identical')

    # First select unique header lines that start with #@.
    PQ_header = []
    for i, header in enumerate(headers):
        for line in header:
            if line.startswith('#@PG'):
                split_line = line.split('\t')
                for j in range(len(split_line)):
                    if (split_line[j].startswith('ID:') 
                        or split_line[j].startswith('PP:')):
                        split_line[j] = split_line[j] + '-' + str(i+1)
                PQ_header.append('\t'.join(split_line))

    other_sam_headers = sum([
        list(set(line for line in header 
            if line.startswith('#@') 
                and (not line.startswith('#@HD'))
                and (not line.startswith('#@SQ'))
                and (not line.startswith('#@PG'))
            ))
        for header in headers], 
        [])
    other_headers = sum([
        list(set(line for line in header 
            if line.startswith('#') 
                and (not line.startswith('#@'))
            ))
        for header in headers], 
        [])

    out_header = [line for line in headers[0] if line.startswith('#@SQ')]
    out_header += PQ_header
    out_header += other_sam_headers
    out_header += other_headers

    out_header = [l.strip() for l in out_header if l.strip()]

    return out_header

if __name__ == '__main__':
    merge()
