import sys
import click

import random

from . import _fileio, _pairsam_format, cli, _headerops, common_io_options

UTIL_NAME = 'pairtools_sample'

@cli.command()

@click.argument(
    'fraction', 
    type=float,
    required=True)

@click.argument(
    'pairs_path', 
    type=str,
    required=False)

@click.option(
    '-o', "--output", 
    type=str, 
    default="", 
    help='output file.'
        ' If the path ends with .gz or .lz4, the output is pbgzip-/lz4c-compressed.'
        ' By default, the output is printed into stdout.')

@click.option(
    '-s', "--seed", 
    type=int, 
    default=None, 
    help='the seed of the random number generator.')

@common_io_options

def sample(
    fraction, pairs_path, output, seed,
    **kwargs
    ):
    '''Select a random subset of pairs in a pairs file.

    FRACTION: the fraction of the randomly selected pairs subset 

    PAIRS_PATH : input .pairs/.pairsam file. If the path ends with .gz or .lz4, the
    input is decompressed by pbgzip/lz4c. By default, the input is read from stdin.

    '''
    sample_py(
        fraction, pairs_path, output, seed,
        **kwargs
    )
    
def sample_py(
    fraction, pairs_path, output, seed,
    **kwargs
    ):

    instream = (_fileio.auto_open(pairs_path, mode='r', 
                                  nproc=kwargs.get('nproc_in'),
                                  command=kwargs.get('cmd_in', None)) 
                if pairs_path else sys.stdin)
    outstream = (_fileio.auto_open(output, mode='w', 
                                   nproc=kwargs.get('nproc_out'),
                                   command=kwargs.get('cmd_out', None)) 
                 if output else sys.stdout)

    header, body_stream = _headerops.get_header(instream)
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    outstream.writelines((l+'\n' for l in header))

    random.seed(seed)

    for line in body_stream:
        if random.random() <= fraction:
            outstream.write(line)

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()


if __name__ == '__main__':
    pair()
