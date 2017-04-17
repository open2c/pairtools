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

CONTEXT_SETTINGS = {
    'help_option_names': ['-h', '--help'],
}

@click.version_option(version=__version__)
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass

from .pairsam_dedup import dedup
from .pairsam_sort import sort
from .pairsam_merge import merge
from .pairsam_markasdup import markasdup
from .pairsam_select import select
from .pairsam_split import split
from .pairsam_restrict import restrict
from .pairsam_parse import parse, parse_cigar, parse_algn
from .pairsam_stats import stats
