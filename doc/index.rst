.. pairtools documentation master file, created by
   sphinx-quickstart on Wed Dec  6 12:32:49 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pairtools
=========

`pairtools` is a simple and fast command-line framework to process sequencing
data from a Hi-C experiment.

`pairtools` process pair-end sequence alignments and perform the following
operations:

- detect ligation junctions (a.k.a. Hi-C pairs) in aligned paired-end sequences of Hi-C DNA molecules
- sort .pairs files for downstream analyses
- detect, tag and remove PCR/optical duplicates 
- generate extensive statistics of Hi-C datasets
- select Hi-C pairs given flexibly defined criteria
- restore .sam alignments from Hi-C pairs

pairtools produce .pairs files compliant with the 
`4DN standards <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_.

The list of available pairtools includes:

.. click_subcommand_table:: pairtools:cli
   :prog: pairtools
   :column_names: Pairtool,Description


Contents:

.. toctree::
   :maxdepth: 3

   quickstart
   installation
   parsing
   sorting
   formats 
   cli_tools


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

