# pairtools

[![Documentation Status](https://readthedocs.org/projects/pairtools/badge/?version=latest)](http://pairtools.readthedocs.org/en/latest/)
[![Build Status](https://travis-ci.org/mirnylab/pairtools.svg?branch=master)](https://travis-ci.org/mirnylab/pairtools)
[![Join the chat at https://gitter.im/mirnylab/distiller](https://badges.gitter.im/mirnylab/distiller.svg)](https://gitter.im/mirnylab/distiller?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1490831.svg)](https://doi.org/10.5281/zenodo.1490831)

## Process Hi-C pairs with pairtools

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

To get started:
- Take a look at a [quick example](https://github.com/mirnylab/pairtools#quick-example)
- Check out the detailed [documentation](http://pairtools.readthedocs.io).

## Data formats

`pairtools` produce and operate on tab-separated files compliant with the
[.pairs](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) 
format defined by the [4D Nucleome Consortium](https://www.4dnucleome.org/). All
pairtools properly manage file headers and keep track of the data
processing history.

Additionally, `pairtools` define the .pairsam format, an extension of .pairs that includes the SAM alignments 
of a sequenced Hi-C molecule. .pairsam complies with the .pairs format, and can be processed by any tool that
operates on .pairs files.

## Installation

Requirements:

- Python 3.x
- Python packages `cython`, `numpy` and `click`.
- Command-line utilities `sort` (the Unix version), `bgzip` (shipped with `tabix`)  and `samtools`. If available, `pairtools` can compress outputs with `pbgzip` and `lz4`.

We highly recommend using the `conda` package manager to install `pairtools` together with all its dependencies. To get it, you can either install the full [Anaconda](https://www.continuum.io/downloads) Python distribution or just the standalone [conda](http://conda.pydata.org/miniconda.html) package manager.

With `conda`, you can install `pairtools` and all of its dependencies from the [bioconda](https://bioconda.github.io/index.html) channel.
```sh
$ conda install -c conda-forge -c bioconda pairtools
```

Alternatively, install `pairtools` and only Python dependencies from PyPI using pip:
```sh
$ pip install pairtools
```

## Quick example

Setup a new test folder and download a small Hi-C dataset mapped to sacCer3 genome:
```bash
$ mkdir /tmp/test-pairtools
$ cd /tmp/test-pairtools
$ wget https://github.com/mirnylab/distiller-test-data/raw/master/bam/MATalpha_R1.bam
```

Additionally, we will need a .chromsizes file, a TAB-separated plain text table describing the names, sizes and the order of chromosomes in the genome assembly used during mapping:
```bash
$ wget https://raw.githubusercontent.com/mirnylab/distiller-test-data/master/genome/sacCer3.reduced.chrom.sizes
```

With `pairtools parse`, we can convert paired-end sequence alignments stored in .sam/.bam format into .pairs, a TAB-separated table of Hi-C ligation junctions:

```bash
$ pairtools parse -c sacCer3.reduced.chrom.sizes -o MATalpha_R1.pairs.gz --drop-sam MATalpha_R1.bam 
```

Inspect the resulting table:

```bash
$ less MATalpha_R1.pairs.gz
```

## Pipelines

- We provide a simple working example of a mapping bash pipeline in /examples/.
- [distiller](https://github.com/mirnylab/distiller-nf) is a powerful
Hi-C data analysis workflow, based on `pairtools` and [nextflow](https://www.nextflow.io/).


## Tools

- `parse`: read .sam files produced by bwa and form Hi-C pairs
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

- `sort`: sort pairs files (the lexicographic order for chromosomes, 
    the numeric order for the positions, the lexicographic order for pair types).

- `merge`: merge sorted .pairs files
    - merge sort .pairs;
    - combine the .pairs headers from all input files;
    - check that each .pairs file was mapped to the same reference genome index 
    (by checking the identity of the @SQ sam header lines).

- `select`: select pairs according to specified criteria
    - select pairs entries according to the provided condition. A programmable
    interface allows for arbitrarily complex queries on specific pair types, 
    chromosomes, positions, strands, read IDs (including matches to a
    wildcard/regexp/list).
    - optionally print the non-matching entries into a separate file.

- `dedup`: remove PCR duplicates from a sorted triu-flipped .pairs file
    - remove PCR duplicates by finding pairs of entries with both sides mapped
    to similar genomic locations (+/- N bp);
    - optionally output the PCR duplicate entries into a separate file.
    - NOTE: in order to remove all PCR duplicates, the input must contain \*all\* 
      mapped read pairs from a single experimental replicate;

- `maskasdup`: mark all pairs in a pairsam as Hi-C duplicates
    - change the field pair_type to DD;
    - change the pair_type tag (Yt:Z:) for all sam alignments;
    - set the PCR duplicate binary flag for all sam alignments (0x400).

- `split`: split a .pairsam file into .pairs and .sam.

- `stats`: calculate various statistics of .pairs files

- `restrict`: identify the span of the restriction fragment forming a Hi-C junction

## Contributing

[Pull requests](https://akrabat.com/the-beginners-guide-to-contributing-to-a-github-project/) are welcome.

For development, clone and install in "editable" (i.e. development) mode with the `-e` option. This way you can also pull changes on the fly.
```sh
$ git clone https://github.com/mirnylab/pairtools.git
$ cd pairtools
$ pip install -e .
```

## License

MIT
