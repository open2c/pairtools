# pairsamtools

[![Documentation Status](https://readthedocs.org/projects/pairsamtools/badge/?version=latest)](http://pairsamtools.readthedocs.org/en/latest/)
[![Build Status](https://travis-ci.org/mirnylab/pairsamtools.svg?branch=master)](https://travis-ci.org/mirnylab/pairsamtools)
[![Join the chat at https://gitter.im/mirnylab/distiller](https://badges.gitter.im/mirnylab/distiller.svg)](https://gitter.im/mirnylab/distiller?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

## build Hi-C mapping pipelines with pairsamtools

pairsamtools is a simple and fast command-line framework to process sequencing
data from a Hi-C experiment.

pairsamtools process pair-end sequence alignments and perform the following
operations:
- detect and classify ligation sites (a.k.a. Hi-C pairs) produced in Hi-C experiments
- sort .pairs files for downstream analyses
- detect, tag and remove PCR/optical duplicates 
- generate extensive statistics of Hi-C datasets
- select Hi-C pairs given flexibly defined criteria
- restore and tag .sam files for selected subsets of Hi-C pairs

To get started, check out the [documentation](http://pairsamtools.readthedocs.io).

pairsamtools produce and operate on tab-separated files compliant with the
[.pairs](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) 
format defined by the [4D Nucleome Consortium](https://www.4dnucleome.org/). All
pairsamtools properly manage file headers and keep track of the data
processing history.

## installation

Requirements:
- python 3.x
- unix sort
- bgzip
- Cython
- numpy
- click

Install using pip:
```sh
$ pip install git+https://github.com/mirnylab/pairsamtools
```

### tools

- parse: read .sam files produced by bwa and form Hi-C pairs
    - form Hi-C pairs by reporting the outer-most mapped positions and the strand
    on the either side of each molecule;
    - report unmapped/multimapped (ambiguous alignments)/chimeric alignments as
    chromosome "!", position 0, strand "-";
    - identify and rescue chrimeric alignments produced by singly-ligated Hi-C 
    molecules with a sequenced ligation junction on one of the sides;
    - perform upper-triangular flipping of the sides of Hi-C molecules 
    such that the first side has a lower sorting index than the second side;
    - form hybrid pairsam output, where each line contains all available data 
    for one Hi-C molecule (outer-most mapped positions on the either side, 
    read ID, pair type, and .sam entries for each alignment);
    - print the .sam header as #-comment lines at the start of the file.

- sort: sort pairsam files (the lexicographic order for chromosomes, 
    the numeric order for the positions, the lexicographic order for pair types).

- merge: merge sorted pairsam files
    - simple merge sort for pairsam entries;
    - combine the pairs headers from all input files;
    - check that each pairsam file was mapped to the same reference genome index 
    (by checking the identity of the @SQ sam header lines).

- select: select pairsam entries with specific field values
    - select pairsam entries according to the provided condition. A programmable
    interface allows for arbitrarily complex queries on specific pair types, 
    chromosomes, positions, strands, read IDs (including matches to a
    wildcard/regexp/list).
    - optionally print the non-matching entries into a separate file.

- dedup: remove PCR duplicates from a sorted triu-flipped pairsam file
    - remove PCR duplicates by finding pairs of entries with both sides mapped
    to similar genomic locations (+/- N bp);
    - optionally output the PCR duplicate entries into a separate file.
    - NOTE: in order to remove all PCR duplicates, the input must contain \*all\* 
      mapped read pairs from a single experimental replicate;

- maskasdup: mark all pairs in a pairsam as Hi-C duplicates
    - change the field pair_type to DD;
    - change the pair_type tag (Yt:Z:) for all sam alignments;
    - set the PCR duplicate binary flag for all sam alignments (0x400).

- split: split a pairsam file into pairs and sam alignments.

- stats: calculate various statistics of .pairs and .pairsam files

- restrict: identify the span of the restriction fragment forming a Hi-C junction

### pipelines

We provide a simple mapping bash pipeline in /examples/.
It serves as an illustration to pairsamtools' functionality and
will not be further developed.

[distiller](https://github.com/mirnylab/distiller-nf) is a powerful
Hi-C data analysis workflow, based on pairsamtools and 
[nextflow](https://www.nextflow.io/).


