Formats for storing Hi-C pairs
==============================

.pairs
------

`.pairs` is a simple tabular format for storing DNA contacts detected in
a Hi-C experiment.  The detailed
`.pairs specification <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_
is defined by the 4DN Consortium.

The body of a .pairs contains a table with a variable number of fields separated by 
a "\\t" character (a horizontal tab). The .pairs specification fixes the content
and the order of the first seven columns:

======== =========== ===============================================================================
 index    name        description  
======== =========== ===============================================================================
 1        read_id     the ID of the read as defined in fastq files 
 2        chrom1      the chromosome of the alignment on side 1 
 3        pos1        the 1-based genomic position of the outer-most (5') mapped bp on side 1 
 4        chrom2      the chromosome of the alignment on side 2 
 5        pos2        the 1-based genomic position of the outer-most (5') mapped bp on side 2 
 6        strand1     the strand of the alignment on side 1 
 7        strand2     the strand of the alignment on side 2 
======== =========== ===============================================================================

A .pairs file starts with a header, an arbitrary number of lines starting
with a "#" character. By convention, the header lines have a format of 
"#field_name: field_value".
The `.pairs specification <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md>`_
mandates a few standard header lines (e.g., column names, 
chromosome order, sorting order, etc), all of which are 
automatically filled in by `pairtools`.

The entries of a .pairs file can be flipped and sorted. "Flipping" means
that *the sides 1 and 2 do not correspond to side1 and side2 in sequencing data.* 
Instead, side1 is defined as the side with the
alignment with a lower sorting index (using the lexographic order for 
chromosome names, followed by the numeric order for positions and the 
lexicographic order for pair types). This particular order of "flipping" is
defined as "upper-triangular flipping", or "triu-flipping". Finally, pairs are
*typically* block-sorted: i.e. first lexicographically by chrom1 and chrom2, 
then numerically by pos1 and pos2.

Pairtools' flavor of .pairs
---------------------------

.pairs files produced by `pairtools` extend .pairs format in a few ways.

1. `pairtools` store null, unmapped, ambiguous (multiply mapped) and chimeric (if not parsed by `parse2` or `--walks-policy all` of `parse`) alignments as chrom='!', pos=0, strand='-'.

#. `pairtools` store the header of the source .sam files in the 
   '#samheader:' fields of the pairs header. When multiple .pairs files are merged,
   the respective '#samheader:' fields are checked for consistency and merged. 

#. Each pairtool applied to .pairs leaves a record in the '#samheader' fields
   (using a @PG sam tag), thus preserving the full history of data processing.

#. `pairtools` append an extra column describing the type of a Hi-C pair:

======== =========== ===============================================================================
 index    name        description  
======== =========== ===============================================================================
 8        pair_type   the type of a Hi-C pair 
======== =========== ===============================================================================

.. _section-pair-types:

Pair types
----------

`pairtools` use a simple two-character notation to define all possible pair
types, according to the quality of alignment of the two sides. The type of a pair 
can be defined unambiguously using the table below. To use this table, 
identify which side has an alignment of a "poorer" quality
(unmapped < multimapped < unique alignment)
and which side has a "better" alignment and find the corresponding row in the table.

======================== ====== =============== ========= ================== ========= ================== ===========
 .                        .      .              Less informative alignment   More informative alignment    .        
------------------------ ------ --------------- ---------------------------- ---------------------------- -----------
Pair type                Code    >2 alignments   Mapped     Unique             Mapped    Unique              Sidedness                           
walk-walk                 WW      |check|         |cross|   |cross|            |cross|   |cross|             0 [1]_
null                      NN      |cross|         |cross|                      |cross|                       0     
corrupt                   XX      |cross|         |cross|                      |cross|                       0 [2]_    
null-multi                NM      |cross|         |cross|                      |check|   |cross|             0     
null-rescued              NR      |check|         |cross|                      |check|   |check|             1 [3]_
null-unique               NU      |cross|         |cross|                      |check|   |check|             1     
multi-multi               MM      |cross|         |check|   |cross|            |check|   |cross|             0     
multi-rescued             MR      |check|         |check|   |cross|            |check|   |check|             1 [3]_
multi-unique              MU      |cross|         |check|   |cross|            |check|   |check|             1     
rescued-unique            RU      |check|         |check|   |check|            |check|   |check|             2 [3]_
unique-rescued            UR      |check|         |check|   |check|            |check|   |check|             2 [3]_
unique-unique             UU      |cross|         |check|   |check|            |check|   |check|             2     
duplicate                 DD      |cross|         |check|   |check|            |check|   |check|             2 [4]_
======================== ====== =============== ========= ================== ========= ================== ===========

.. [1] "walks", or, `C-walks <https://www.nature.com/articles/nature20158>`_ are
   Hi-C molecules formed via multiple ligation events which cannot be reported 
   as a single pair.  

.. [2] "corrupt" pairs are those with technical issues - e.g. missing a 
   FASTQ sequence/SAM entry from one side of the molecule.

.. [2] "rescued" pairs have two non-overlapping alignments on one of the sides
   (referred below as the chimeric side/read), but the inner (3'-) one extends the 
   only alignment on the other side (referred as the non-chimeric side/read).
   Such pairs form when one of the two ligated DNA fragments is shorter than
   the read length. In this case, one of the reads contains this short fragment
   entirely, together with the ligation junction and a chunk of the other DNA fragment 
   (thus, this read ends up having two non-overlapping alignments).
   Following the procedure introduced in `HiC-Pro <https://github.com/nservant/HiC-Pro>`_
   and `Juicer <https://github.com/theaidenlab/juicer>`_, `pairtools parse` 
   rescues such Hi-C molecules, reports the position of the 5' alignment on the
   chimeric side, and tags them as "NU", "MU", "UR" or "RU" pair type, depending 
   on the type of the 5' alignment on the chimeric side. Such molecules can and
   should be used in downstream analysis.
   Read more on the rescue procedure in :doc:`the section on parsing <parsing>`.

.. [3] `pairtools dedup` detects molecules that could be formed via PCR duplication and
   tags them as "DD" pair type. These pairs should be excluded from downstream 
   analyses.

.pairsam 
--------

`pairtools` also define .pairsam, a valid extension of the .pairs format.
On top of the pairtools' flavor of .pairs, .pairsam format adds two extra 
columns containing the alignments from which the Hi-C pair was extracted:

======== =========== ===============================================================================
 index    name        description  
======== =========== ===============================================================================
 9        sam1        the sam alignment(s) on side 1; separate supplemental alignments by NEXT_SAM
 10       sam2        the sam alignment(s) on side 2; separate supplemental alignments by NEXT_SAM
======== =========== ===============================================================================

Note that, normally, the fields of a sam alignment are separated by a horizontal 
tab character (\\t), which we already use to separate .pairs columns. To
avoid confusion, we replace the tab character in sam entries stored in sam1 and 
sam2 columns with a UNIT SEPARATOR character (\\031).

Finally, sam1 and sam2 can store multiple .sam alignments, separated by a string
'\\031NEXT_SAM\\031'


.. |check| unicode:: U+2714 .. check
.. |cross| unicode:: U+274C .. cross

Extra columns
----------------

`pairtools` can operate on `.pairs/.pairsam` with extra columns.
Extra columns are specified in the order defined by the order their addition by various tools.
Column names can be checked in the header of `.pairs/.pairsam` file.
We provide `pairtools header` utilities for manipulating and verifying compatibility of headers and their columns.

The list of additional columns used throughout `pairtools` modules:

=================================== =================== ====================== ================================================== =================
extra column                        generating module   format                 how to add it                                       description
=================================== =================== ====================== ================================================== =================
mapq1, mapq2                         `parse/parse2`      number from 0 to 255   `pairtools parse --add-columns mapq`              `Mapping quality <https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess>`_, as reported in .sam/.bam, $-10 log_{10}(P_{error})$
pos51, pos52                         `parse/parse2`      genomic coordinate     `pairtools parse --add-columns pos5`              5' position of alignment (closer to read start)
pos31, pos32                         `parse/parse2`      genomic coordinate     `pairtools parse --add-columns pos3`              3' position of alignment (further from read start)
cigar1, cigar2                       `parse/parse2`      string                 `pairtools parse --add-columns cigar`             `CIGAR, or Compact Idiosyncratic Gapped Alignment Report <https://en.wikipedia.org/wiki/Sequence_alignment#CIGAR_Format>`_ of alignment, as reported in .sam/.bam
read_len1, read_len2                 `parse/parse2`      number                 `pairtools parse --add-columns read_len`          read length
matched_bp1, matched_bp2             `parse/parse2`      number                 `pairtools parse --add-columns matched_bp`        number of matched alignment basepairs to the reference
algn_ref_span1, algn_ref_span2       `parse/parse2`      number                 `pairtools parse --add-columns algn_ref_span`     basepairs of reference covered by alignment
algn_read_span1, algn_read_span2     `parse/parse2`      number                 `pairtools parse --add-columns algn_read_span`    basepairs of read covered by alignment
dist_to_51, dist_to_52               `parse/parse2`      number                 `pairtools parse --add-columns dist_to_5`         distance to 5'-end of read
dist_to_31, dist_to_32               `parse/parse2`      number                 `pairtools parse --add-columns dist_to_3`         distance to 3'-end of read
seq1, seq2                           `parse/parse2`      string                 `pairtools parse --add-columns seq`               sequence of alignment
mismatches1, mismatches2             `parse/parse2`      string                 `pairtools parse --add-columns mismatches`        comma-separated list of mismatches relative to the reference, "{ref_letter}:{mut_letter}:{phred}:{ref_position}:{read_position}"
XB1/2,AS1/2,XS1/2 or any sam tag     `parse/parse2`                             `pairtools parse --add-columns XA,XB,NM`          format depends on `tag specification <https://samtools.github.io/hts-specs/SAMv1.pdf>`_
walk_pair_type                       `parse/parse2`      string                 `pairtools parse2 --add-pair-index`               Type of the pair relative to R1 and R2 reads of paired-end sequencing, see `pasring docs <https://pairtools.readthedocs.io/en/latest/parsing.html#rescuing-complex-walks>`_
walk_pair_index                      `parse/parse2`      number                 `pairtools parse2 --add-pair-index`               Order of the pair in the complex walk, starting from 5'-end of left read, see `pasring docs <https://pairtools.readthedocs.io/en/latest/parsing.html#rescuing-complex-walks>`_
phase                                `phase`             0, 1 or "."            `pairtools phase`                                 Phase of alignment (haplotype 1, 2, on unphased), see `phasing walkthrough <https://pairtools.readthedocs.io/en/latest/examples/pairtools_phase_walkthrough.html>`_
rfrag1, rfrag2                       `restrict`          number                 `pairtools restrict`                              Unique index of the restriction fragment after annotating pairs positions, see `restriction walkthrough <https://pairtools.readthedocs.io/en/latest/examples/pairtools_restrict_walkthrough.html>`_
rfrag_start1, rfrag_start2           `restrict`          number                 `pairtools restrict`                              Coordinate of the start of restriction fragment
rfrag_end1, rfrag_end2               `restrict`          number                 `pairtools restrict`                              Coordinate of the end of restriction fragment
=================================== =================== ====================== ================================================== =================

