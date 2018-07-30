.. pairtools documentation master file, created by
   sphinx-quickstart on Wed Dec  6 12:32:49 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Overview
========

`pairtools` is a simple and fast command-line framework to process sequencing
data from a Hi-C experiment. `pairtools` perform various operations on Hi-C
pairs and occupy the middle position in a typical Hi-C data processing 
pipeline:

.. figure:: _static/hic-processing-pipeline.png
   :width: 100%
   :alt: The diagram of a typical processing pipeline for Hi-C data
   :align: center

   In a typical Hi-C pipeline, DNA sequences (reads) are aligned to the 
   reference genome, converted into ligation junctions and binned, thus 
   producing a Hi-C contact map.


`pairtools` aim to be an all-in-one tool for processing Hi-C pairs, and 
can perform following operations:

- detect ligation junctions (a.k.a. Hi-C pairs) in aligned paired-end sequences of Hi-C DNA molecules
- sort .pairs files for downstream analyses
- detect, tag and remove PCR/optical duplicates 
- generate extensive statistics of Hi-C datasets
- select Hi-C pairs given flexibly defined criteria
- restore .sam alignments from Hi-C pairs

`pairtools` produce .pairs files compliant with the 
`4DN standard <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_.

The full list of available pairtools:

.. click_subcommand_table:: pairtools:cli
   :prog: pairtools
   :column_names: Pairtool,Description


Contents:

.. toctree::
   :hidden:

   self

.. toctree::
   :maxdepth: 3

   quickstart
   installation
   parsing
   sorting
   formats 
   technotes
   cli_tools


* :ref:`genindex`

