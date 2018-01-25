# -*- coding: utf-8 -*-
"""
pairsamtools
~~~~~~~~~~~~

CLI tools to process mapped Hi-C data

:copyright: (c) 2017 Massachusetts Institute of Technology
:author: Mirny Lab
:license: MIT

"""

__version__ = '0.0.1-dev'


import click
import functools
import sys

CONTEXT_SETTINGS = {
    'help_option_names': ['-h', '--help'],
}


@click.version_option(version=__version__)
@click.group(context_settings=CONTEXT_SETTINGS)
@click.option(
    '--post-mortem', 
    help="Post mortem debugging", 
    is_flag=True,
    default=False
)

@click.option(
    '--output-profile', 
    help="Profile performance with Python cProfile and dump the statistics "
         "into a binary file", 
    type=str,
    default=''
)
def cli(post_mortem, output_profile):
    if post_mortem:
        import traceback
        try:
            import ipdb as pdb
        except ImportError:
            import pdb
        def _excepthook(exc_type, value, tb):
            traceback.print_exception(exc_type, value, tb)
            print()
            pdb.pm()
        sys.excepthook = _excepthook

    if output_profile:
        import cProfile 
        import atexit
        
        pr = cProfile.Profile()
        pr.enable()

        def _atexit_profile_hook():
            pr.disable()
            pr.dump_stats(output_profile)

        atexit.register(_atexit_profile_hook)


def common_io_options(func):
    @click.option(
        '--nproc-in',
        type=int, 
        default=3, 
        show_default=True,
        help='Number of processes used by the auto-guessed input decompressing command.'
        )
    @click.option(
        '--nproc-out',
        type=int, 
        default=8, 
        show_default=True,
        help='Number of processes used by the auto-guessed output compressing command.'
        )
    @click.option(
        '--cmd-in',
        type=str, 
        default=None, 
        help='A command to decompress the input file. '
             'If provided, fully overrides the auto-guessed command. '
             'Does not work with stdin. '
             'Must read input from stdin and print output into stdout. '
             'EXAMPLE: pbgzip -dc -n 3'
        )
    @click.option(
        '--cmd-out',
        type=str, 
        default=None, 
        help='A command to compress the output file. '
             'If provided, fully overrides the auto-guessed command. '
             'Does not work with stdout. '
             'Must read input from stdin and print output into stdout. '
             'EXAMPLE: pbgzip -c -n 8'
        )

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper

from .pairsam_dedup import dedup
from .pairsam_sort import sort
from .pairsam_merge import merge
from .pairsam_markasdup import markasdup
from .pairsam_select import select
from .pairsam_split import split
from .pairsam_restrict import restrict
from .pairsam_parse import parse, parse_cigar, parse_algn
from .pairsam_stats import stats
