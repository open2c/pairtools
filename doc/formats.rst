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

1. `pairtools` store null/ambiguous/chimeric alignments as chrom='!', pos=0, strand='-'.

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

Pair types
----------

`pairtools` use a simple two-character notation to define all possible pair
types, according to the quality of alignment of the two sides. The type of a pair 
can be defined unambiguously using the table below. To use this table, 
identify which side has an alignment of a "poorer" quality
(unmapped < multimapped < unique alignment)
and which side has a "better" alignment and find the corresponding row in the table.

=============== ========= ==================== ========= =================== ======================== ========== ===========
  .              Less informative alignment     More informative alignment    .                        .          .        
--------------- ------------------------------ ----------------------------- ------------------------ ---------- -----------
 >2 alignments   Mapped    Unique               Mapped    Unique              Pair type                Code       Sidedness                           
 |check|         |cross|   |cross|              |cross|   |cross|             chimeric-chimeric        CC         0 [1]_
 |cross|         |cross|                        |cross|                       null                     NN         0     
 |cross|         |cross|                        |check|   |cross|             null-multi               NM         0     
 |cross|         |cross|                        |check|   |check|             null-unique              NU         1     
 |check|         |cross|                        |check|   |check|             null-rescued-chimeric    NR         1 [2]_
 |cross|         |check|   |cross|              |check|   |cross|             multi-multi              MM         0     
 |cross|         |check|   |cross|              |check|   |check|             multi-unique             MU         1     
 |check|         |check|   |cross|              |check|   |check|             multi-rescued-chimeric   MR         2 [2]_
 |cross|         |check|   |check|              |check|   |check|             unique-unique            UU         2     
 |check|         |check|   |check|              |check|   |check|             rescured-chimeric        UR or RU   2 [2]_
 |cross|         |check|   |check|              |check|   |check|             duplicate                DD         2 [3]_
=============== ========= ==================== ========= =================== ======================== ========== ===========

.. [1] chimeric reads represent Hi-C molecules formed via multiple ligation
   events and thus cannot be reported as a single pair.

.. [2] some chimeric reads correspond to valid Hi-C molecules formed via a single
   ligation event, with the ligation junction sequenced through on one side. 
   Following the procedure introduced in `HiC-Pro <https://github.com/nservant/HiC-Pro>`_
   and `Juicer <https://github.com/theaidenlab/juicer>`_, `pairtools` rescue such 
   molecules, report their outer-most mapped positions and tag them as "UR" or "RU" pair type.
   Such molecules can and should be used in downstream analysis.

.. [3] pairtools detect molecules that could be formed via PCR duplication and
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

Technical notes
---------------

The motivation behind some of the technical decisions in the definition
of .pairsam:

- while the information in columns 1-8 may appear redundant to sam alignments in
  the columns 9+, extracting this information is non-trivial and thus is better 
  done only once with results stored.
- storing sam entries together with pairs drastically speeds up and simplifies 
  several operations like filtering and tagging of unmapped/ambiguous/duplicated 
  Hi-C molecules.
- the exclamation mark "!" is used as a character for unmapped chromosomes
  because it has a lexicographic sorting order lower than that of "0", good 
  interpretability and no other reserved technical roles.



.. |check| unicode:: U+2714 .. check
.. |cross| unicode:: U+274C .. cross

