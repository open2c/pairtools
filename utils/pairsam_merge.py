#!/usr/bin/env python
import sys
import glob
import subprocess
import argparse
from _distiller_common import open_bgzip, DISTILLER_VERSION, \
    append_pg_to_sam_header, get_header


def main():
    parser = argparse.ArgumentParser(
        description="""Merge multiple sorted pairsam files. 
        The @SQ records of the SAM header must be identical; the sorting order of 
        these lines is taken from the first file in the list. 
        The ID fields of the @PG records of the SAM header are modified with a
        numeric suffix to produce unique records.
        The other unique SAM and non-SAM header lines are copied into the output header.
        """
    )
    parser.add_argument(
        'infile', 
        nargs='+', 
        type=str,
        help='a file to merge or a group of files specified by a wildcard',
        default=sys.stdin)

    args = vars(parser.parse_args())
    paths = sum([glob.glob(mask) for mask in args['infile']], [])
    merged_header = form_merged_header(paths)

    merged_header = append_pg_to_sam_header(
        merged_header,
        {'ID': 'pairsam_merge',
         'PN': 'pairsam_merge',
         'VN': DISTILLER_VERSION,
         'CL': ' '.join(sys.argv)
         })

    for line in merged_header:
        print(line, end='', flush=True)

    command = r'''/bin/bash -c 'sort -k 1,1 -k 4,4 -k 2,2n -k 5,5n -k 8,8 --field-separator=$'\''\v'\'' --merge'''
    for path in paths:
        if path.endswith('.gz'):
            command += r''' <(zcat {} | sed -n -e '\''/^[^#]/,$p'\'')'''.format(path)
        else:
            command += r''' <(sed -n -e '\''/^[^#]/,$p'\'' {})'''.format(path)
    command += "'"
    subprocess.call(command, shell=True)


def form_merged_header(paths):
    headers = []
    for path in paths:
        header = []
        f = open_bgzip(path, mode='r')
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
    main()
