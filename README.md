# distiller

## map Hi-C reads with distiller

distiller is a simple and fast command-line framework to process sequencing
information from a Hi-C experiment.

distiller provides command-line tools to perform the following operations:
- map the sequences of Hi-C molecules to a reference genome
- form and classify Hi-C pairs, accounting for the ligation junctions
- filter mapped sequences for valid Hi-C pairs
- remove PCR duplicates 
- tag .sam entries with the Hi-C specific information
- produce sorted lists of Hi-C pairs for downstream analyses

distiller provides a simple bash pipeline.

## installation

At the moment, distiller is provided as a set of command-line scripts.
In the nearest future, distiller will become pip/conda-installable.

Requirements:
- python 3.x
- bgzip
- Cython
- bwa

## usage

### pipeline



### tools

- sam_to_pairsam: read .sam files produced by bwa and form Hi-C pairs
    - 
    - form pairsam output with

## data conventions

### pairsam

### pair types


| read1  |                 |                                 | read2  |                 |                                 | Code     | Pair type         | Sidedness |
|--------|-----------------|---------------------------------|--------|-----------------|---------------------------------|----------|-------------------|-----------|
| Mapped | Uniquely mapped | Linear (non-chimeric) alignment | Mapped | Uniquely mapped | Linear (non-chimeric) alignment |          |                   |           |
| -      |                 |                                 | -      |                 |                                 | NN       | null              | 0         |
| -      |                 |                                 | +      | -               |                                 | NM       | null-multi        | 0         |
| -      |                 |                                 | +      | +               | -                               | NC       | null-chimeric     | 0*        |
| -      |                 |                                 | +      | +               | +                               | NL       | null-linear       | 1         |
| +      | -               |                                 | +      | -               |                                 | MM       | multi-multi       | 0         |
| +      | -               |                                 | +      | +               | -                               | MC       | multi-chimeric    | 0*        |
| +      | -               |                                 | +      | +               | +                               | ML       | multi-linear      | 1         |
| +      | +               | -                               | +      | +               | -                               | CC       | chimeric-chimeric | 0*        |
| +      | +               | -                               | +      | +               | +                               | CL or CX | chimeric-linear   | 1* or 2**   |
| +      | +               | +                               | +      | +               | +                               | LL       | linear-linear     | 2         |
| +      | +               | +                               | +      | +               | +                               | DD       | duplicate         | 2***         |

* - chimeric reads may represent Hi-C molecules formed via multiple ligation
events and thus cannot be interpreted as unambigous pair

** - some chimeric reads correspond to valid Hi-C molecules formed via a single
ligation event, with the ligation junction sequenced through on one side. 
Following the procedure introduced in
[Juicer](https://github.com/theaidenlab/juicer), distiller rescues such 
molecules, reports their 5' mapped positions and tags them as "CX" pair type.
Such molecules can and should be used in downstream analysis.

*** - distiller detects molecules that could be formed via PCR duplication and
tags them as "DD" pair type. These pairs should be excluded from downstream 
analyses.
