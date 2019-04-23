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

=============== ========= ================== ========= ================== ======================== ====== ===========
  .              Less informative alignment   More informative alignment   .                        .      .        
--------------- ---------------------------- ---------------------------- ------------------------ ------ -----------
 >2 alignments   Mapped    Unique             Mapped    Unique             Pair type                Code   Sidedness                           
 |check|         |cross|   |cross|            |cross|   |cross|            walk-walk                WW     0 [1]_
 |cross|         |cross|                      |cross|                      null                     NN     0     
 |cross|         |cross|                      |cross|                      corrupt                  XX     0 [2]_    
 |cross|         |cross|                      |check|   |cross|            null-multi               NM     0     
 |check|         |cross|                      |check|   |check|            null-rescued             NR     1 [3]_
 |cross|         |cross|                      |check|   |check|            null-unique              NU     1     
 |cross|         |check|   |cross|            |check|   |cross|            multi-multi              MM     0     
 |check|         |check|   |cross|            |check|   |check|            multi-rescued            MR     1 [3]_
 |cross|         |check|   |cross|            |check|   |check|            multi-unique             MU     1     
 |check|         |check|   |check|            |check|   |check|            rescued-unique           RU     2 [3]_
 |check|         |check|   |check|            |check|   |check|            unique-rescued           UR     2 [3]_
 |cross|         |check|   |check|            |check|   |check|            unique-unique            UU     2     
 |cross|         |check|   |check|            |check|   |check|            duplicate                DD     2 [4]_
=============== ========= ================== ========= ================== ======================== ====== ===========

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

