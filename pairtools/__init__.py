# -*- coding: utf-8 -*-
"""
pairtools
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

from .dedup import dedup
from .sort import sort
from .merge import merge
from .markasdup import markasdup
from .select import select
from .extractsam import extractsam
from .restrict import restrict
from .parse import parse, parse_cigar, parse_algn
from .stats import stats
