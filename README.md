# distiller

[![Build Status](https://travis-ci.org/mirnylab/distiller.svg?branch=master)](https://travis-ci.org/mirnylab/distiller)
[![Join the chat at https://gitter.im/mirnylab/distiller](https://badges.gitter.im/mirnylab/distiller.svg)](https://gitter.im/mirnylab/distiller?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

## map Hi-C reads with distiller

distiller is a simple and fast command-line framework to process sequencing
data from a Hi-C experiment.

distiller provides command-line tools to perform the following operations:
- map the sequences of Hi-C molecules to a reference genome
- form and classify Hi-C pairs, accounting for the ligation junctions
- filter valid Hi-C pairs
- remove PCR duplicates 
- tag .sam entries with Hi-C specific information
- produce sorted lists of Hi-C pairs for downstream analyses

distiller is provided as a set of command-line tools and a simple example bash
pipeline. In the nearest future, distiller will include a simple make-like 
workflow.

## installation

At the moment, distiller is provided as a set of command-line scripts.
In the nearest future, distiller will become pip/conda-installable.

Requirements:
- python 3.x
- bgzip
- bwa
- Cython
- numpy
- click

## usage

### tools

- sam_to_pairsam: read .sam files produced by bwa and form Hi-C pairs
    - form Hi-C pairs by reporting the outer-most mapped positions and the strand
    on the either side of each molecule;
    - report unmapped/multimapped (ambiguous alignments)/chimeric alignments as
    chromosome "!", position 0, strand "-";
    - find and rescue chrimeric alignments produced by singly-ligated Hi-C 
    molecules with a sequenced ligation junction on one of the sides;
    - perform upper-triangular flipping of the sides of Hi-C molecules 
    such that the first side has a lower sorting index than the second side;
    - form hybrid pairsam output, where each line contains all available data 
    for one Hi-C molecule (outer-most mapped positions on the either side, 
    read ID, pair type, and .sam entries for each alignment);
    - print the .sam header as #-comment lines at the start of the file.

- pairsam_sort: sort pairsam files (the lexicographic order for chromosomes, 
    the numeric order for the positions, the lexicographic order for pair types).

- pairsam_merge: merge sorted pairsam files
    - simple merge sort for pairsam entries;
    - combine the #-comment sections from the beginning of each file. Report the
    sam header @SQ lines first, then other sam header lines, then non-sam
    comments;
    - check that each pairsam file was mapped to the same reference genome index 
    (by checking the identity of the @SQ sam header lines).

- pairsam_select: select pairsam entries with specific field values
    - select pairsam entries with specific pair types, chromosomes or
    read IDs (allow matching to a wildcard/regexp/list).
    - optionally print the non-matching entries into a separate file.

- pairsam_dedup: remove PCR duplicates from a sorted triu-flipped pairsam file
    - remove PCR duplicates by finding pairs of entries with both sides mapped
    to similar genomic locations (+/- N bp);
    - optionally output the PCR duplicate entries into a separate file.
    - NOTE: in order to remove all PCR duplicates, the input must contain \*all\* 
    LL/CX read pairs from a single experimental replicate;

- pairsam_maskasdup: mark all pairs in a pairsam as Hi-C duplicates
    - change the field pair_type to DD;
    - change the pair_type tag (Yt:Z:) for all sam alignments;
    - set the PCR duplicate binary flag for all sam alignments (0x400).

- pairsam_split: split a pairsam file into pairs and sam alignments.

- pairsam_get_header: print the header of a pairsam file

- pairsam_skip_header: print the body of a pairsam file, skipping the header

### pipeline

Currently, distiller provides a simple mapping bash pipeline in /examples/.
This bash pipeline serves as an illustration to distiller's functionality and
will not be further developed.

In the nearest future, distiller will receive a make-like workflow for flexible
and reliable execution.

## data conventions

### pairsam

distiller defines pairsam, a simple tabular format to pass and process
alignments of Hi-C molecules.

A pairsam starts with an arbitrary number of header lines, each starting with
a "#" character.

The body of a pairsam contains a table with a variable number of fields separated by 
a "\v" character (a vertical tab):

| index | name      | description |
|-------|-----------|-------------|
| 1     | read_id   | the ID of the read as defined in fastq files |
| 2     | chrom1    | the chromosome of the alignment on side 1 |
| 3     | chrom2    | the chromosome of the alignment on side 2 |
| 4     | pos1      | the 1-based genomic position of the outer-most (5') mapped bp on side 1 |
| 5     | pos2      | the 1-based genomic position of the outer-most (5') mapped bp on side 2 |
| 6     | strand1   | the strand of the alignment on side 1 |
| 7     | strand2   | the strand of the alignment on side 2 |
| 8     | pair_type | the type of a Hi-C pair |
| 9     | sam1      | the sam alignment(s) on side 1; separate supplemental alignments by NEXT_SAM|
| 10    | sam2      | the sam alignment(s) on side 2; separate supplemental alignments by NEXT_SAM|

*The sides 1 and 2 as defined in pairsam file do not correspond to side1 and
side2 in sequencing data!* Instead, side1 is defined as the side with the
alignment with a lower sorting index (using the lexographic order for 
chromosome names, followed by the numeric order for positions and the 
lexicographic order for pair types). This procedure is defined as 
upper-triangular flipping, or triu-flipping.

The rows of the table are block-sorted: i.e. first lexicographically 
by chrom1 and chrom2, then numerically by pos1 and pos2, then lexicographically
by pair_type.

Null/ambiguous/chimeric alignments are stored as chrom='!', pos=0, strand='-'.

Notes of the motivation behind some of the technical decisions in the definition
of pairsam:
- while the information in columns 1-8 may appear redundant to sam alignments in
the columns 9+, extracting this information is non-trivial and thus is better 
done only once with results stored.
- storing sam entries together with pairs drastically speeds up and simplifies 
several operations like filtering and tagging of unmapped/ambiguous/duplicated 
Hi-C molecules.
- pair flipping and sorting is essential for the processing steps like PCR
duplicate removal and aggregation.
- the vertical tab is used as a field separator because the horizontal tab is 
already used by sam entries. All printable characters can apprear in the
PHRED field of sam entries and thus cannot be used as a separator; out of all 
control characters, the vertical tab is the only remaining unused character 
with a common escape notation.
- the exclamation mark "!" is used as a character for unmapped chromosomes
because it has a lexicographic sorting order lower than that of "0", good 
interpretability and no other reserved technical roles.

### pair types

distiller uses a simple two-character notation to define all possible pair types
by the quality of alignment. For each pair, its type can be defined unambigously
using the table below. To use this table, indentify which side has an alignment 
of a "poorer" quality (unmapped < multimapped < chimeric alignment < linear alignment)
and which side has a "better" alignment and find the corresponding row in the table.




| more poorly aligned side  |                 |                                 | better aligned side  |                 |                                 | Code     | Pair type         | Sidedness |
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

\*  chimeric reads may represent Hi-C molecules formed via multiple ligation
events and thus cannot be interpreted as unambigous pairs.

**  some chimeric reads correspond to valid Hi-C molecules formed via a single
ligation event, with the ligation junction sequenced through on one side. 
Following the procedure introduced in
[Juicer](https://github.com/theaidenlab/juicer), distiller rescues such 
molecules, reports their outer-most mapped positions and tags them as "CX" pair type.
Such molecules can and should be used in downstream analysis.

***  distiller detects molecules that could be formed via PCR duplication and
tags them as "DD" pair type. These pairs should be excluded from downstream 
analyses.
