.. pairtools documentation master file, created by
   sphinx-quickstart on Wed Dec  6 12:32:49 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pairtools
=========

pairtools is a set of simple and fast command-line tools to process 
sequencing data from Hi-C experiments.

pairtools operate on sequence alignments and perform the following operations:

* detect and classify ligation sites (a.k.a. `Hi-C pairs`) produced in Hi-C experiments
* sort Hi-C pairs for downstream analyses 
* detect, tag and remove PCR/optical duplicates
* generate extensive statistics of Hi-C datasets
* select Hi-C pairs given flexibly defined criteria
* restore and tag .sam files for selected subsets of Hi-C pairs

pairtools produce .pairs files compliant with the 
`4DN standards <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_.


Contents:

.. toctree::
   :maxdepth: 3

   quickstart
   parsing
   pairsam
   cli_tools


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

