Parsing sequence alignments into Hi-C pairs
===========================================

Overview
--------

Hi-C experiments aim to measure the frequencies of contacts between all pairs
of loci in the genome. In these experiments, the spacial structure of chromosomes 
is first fixed with formaldehyde crosslinks, after which DNA is partially
digested with restriction enzymes and then re-ligated back. Then, DNA is 
shredded into smaller pieces, released from the nucleus, sequenced and aligned to
the reference genome. The resulting sequence alignments reveal if DNA molecules 
were formed through ligations between DNA from different locations in the genome.
These ligation events imply that ligated loci were close to each other
when the ligation enzyme was active, i.e. they formed "a contact".

``pairtools parse`` detects ligation events in the aligned sequences of 
DNA molecules formed in Hi-C experiments and reports them in the .pairs/.pairsam 
format.

Terminology 
-----------

Throughout this document we will be using the same visual language to describe
how DNA sequences (in the .fastq format) are transformed into sequence alignments 
(.sam/.bam) and into ligation events (.pairs).

.. figure:: _static/terminology.png
   :scale: 50 %
   :alt: The visual language to describe transformation of Hi-C data
   :align: center

   DNA sequences (reads) are aligned to the reference genome and converted into
   ligation events

Short-read sequencing determines the sequences of the both ends (or, **sides**)
of DNA molecules (typically 50-300 bp), producing **read pairs** in .fastq format 
(shown in the first row on the figure above).
In such reads, base pairs are reported from the tips inwards, which is also
defined as the **5'->3'** direction (in accordance of the 5'->3' direction of the
DNA strand that sequence of the corresponding side of the read).

Alignment software maps both reads of a pair to the reference genome, producing
**alignments**, i.e. segments of the reference genome with matching sequences.
Typically, if the read length is not very large (< 150 bp), there will be only
two alignments per read pair, one on each side. But, sometimes, the parts of one
or both sides may map to different locations on the genome, producing more than
two alignments per DNA molecule (see :ref:`section-walks`).

``pairtools parse`` converts alignments into **ligation events** (aka
**Hi-C pairs** aka **pairs**). In the simplest case, when each side has only one 
unique alignment (i.e. the whole side maps to a single unique segment of the 
genome), for each side, we report the chromosome, the genomic position of the
outer-most (5') aligned base pair and the strand of the reference genome that 
the read aligns to.  ``pairtools parse`` assigns to such pairs the type ``UU``
(unique-unique).

Unmapped/multimapped reads
--------------------------

Sometimes, one side or both sides of a read pair may not align to the 
reference genome:

.. figure:: _static/read_pair_NU_NN.png
   :scale: 50 %
   :alt: Read pairs missing an alignment on one or both sides
   :align: center

   Read pairs missing an alignment on one or both sides

In this case, ``pairtools parse`` fills in the chromosome of the corresponding
side of Hi-C pair with ``!``, the position with ``0`` and the strand with ``-``.
Such pairs are reported as type ``NU`` (null-unique, when the other side has
a unique alignment) or ``NN`` (null-null, when both sides lack any alignment).

Similarly, when one or both sides map to many genome locations equally well (i.e.
have non-unique, or, multi-mapping alignments), ``pairtools parse`` reports 
the corresponding sides as (chromosome= ``!``, position= ``0``, strand= ``-``) and 
type ``MU`` (multi-unique) or ``MM`` (multi-multi) or ``NM`` (null-multi),
depending on the type of the alignment on the other side.

.. figure:: _static/read_pair_MU_MM_NM.png
   :scale: 50 %
   :alt: Read pairs with a non-unique alignment on one or both sides
   :align: center

   Read pairs with a non-unique (multi-) alignment on one side
   
``pairtools parse`` calls an alignment to be multi-mapping when its
`MAPQ score <https://bioinformatics.stackexchange.com/questions/2417/meaning-of-bwa-mem-mapq-scores>`_
(which depends on the scoring gap between the two best candidate alignments for a segment)
is equal or greater than the value specied with the ``--min-mapq`` flag (by default, 1).

.. _section-walks:

Multiple ligations (walks)
--------------------------

If the read is long enough (e.g. larger than 150 bp), it may contain more than two alignments:

.. figure:: _static/read_pair_WW.png
   :scale: 50 %
   :alt: A sequenced Hi-C molecule that was formed via multiple ligations
   :align: center

   A sequenced Hi-C molecule that was formed via multiple ligations

Molecules like these typically form via multiple ligation events and we call them
walks [1]_. The mode of walks reporting is controlled by ``--walks-policy`` parameter of ``pairtools parse``.
You can report all the alignments in the reads by using ``pairtools parse2`` (see :ref:`parse2`).

A pair of sequential alignments on a single read is **ligation junction**. Ligation junctions are the Hi-C contacts
that have been directly observed in the experiment. However, traditional Hi-C pairs do not have direct evidence of ligation
because they arise from read pairs that do not necessarily contain ligation junction.
To filter out the molecules with complex walks, ``--walks-policy`` can be set to:

- ``mask`` to tag these molecules as type ``WW`` (single ligations are rescued, see :ref:`section-single-ligation-rescue`) ,
- ``5any`` to report the 5'-most alignment on each side,
- ``5unique`` to report the 5'-most unique alignment on each side,
- ``3any`` to report the 3'-most alignment on each side,
- ``3unique`` to report the 3'-most unique alignment on each side.

.. _section-complex-walks-rescue:

Rescuing single ligations
-------------------------

Importantly, some of DNA molecules containing only one ligation junction
may still end up with three alignments:

.. figure:: _static/read_pair_UR.png
   :scale: 50 %
   :alt: Not all read pairs with three alignments come from "walks"
   :align: center

   Not all read pairs with three alignments come from "walks"

A molecule formed via a single ligation gets three alignments when one of the 
two ligated DNA pieces is shorter than the read length, such that that read on 
the corresponding side sequences through the ligation junction and into the other 
piece [2]_. The amount of such molecules depends on the type of the restriction 
enzyme, the typical size of DNA molecules in the Hi-C library and the read 
length, and sometimes can be considerable.

``pairtools parse`` detects such molecules and **rescues** them (i.e.
changes their type from a *walk* to a single-ligation molecule). It tests
walks with three aligments using three criteria:

.. figure:: _static/read_pair_UR_criteria.png
   :scale: 50 %
   :alt: The three criteria used for "rescue"
   :align: center

   The three criteria used to "rescue" three-alignment walks: cis, point towards each other, short distance

1. On the side with two alignments (the **chimeric** side), the "inner" (or, 3') 
   alignment must be on the same chromosome as the alignment on the non-chimeric
   side.

2. The "inner" alignment on the chimeric side and the alignment on the 
   non-chimeric side must point toward each other.

3. These two alignments must be within the distance specified with the
   ``--max-molecule-size`` flag (by default, 2000bp).

Sometimes, the "inner" alignment on the chimeric side can be non-unique or "null" 
(i.e. when the unmapped segment is longer than ``--max-inter-align-gap``, 
as described in :ref:`section-gaps`). ``pairtools parse`` ignores such alignments
altogether and thus rescues such *walks* as well.

.. figure:: _static/read_pair_UR_MorN.png
   :scale: 50 %
   :alt: A walk with three alignments get rescued, when the middle alignment is multi- or null
   :align: center

   A walk with three alignments get rescued, when the middle alignment is multi- or null.


.. _section-gaps:

Interpreting gaps between alignments
------------------------------------

Reads that are only partially aligned to the genome can be interpreted in
two different ways. One possibility is to assume that this molecule
was formed via at least two ligations (i.e. it's a *walk*) but the non-aligned
part (a **gap**) was missing from the reference genome for one reason or another.
Another possibility is to simply ignore this gap (for example, because it could
be an insertion or a technical artifact), thus assuming that our
molecule was formed via a single ligation and has to be reported:

.. figure:: _static/read_pair_gaps_vs_null_alignment.png
   :scale: 50 %
   :alt: A gap between alignments can be ignored or interpeted as a "null" alignment
   :align: center

   A gap between alignments can interpeted as a legitimate segment without
   an alignment or simply ignored

Both options have their merits, depending on a dataset, quality of the reference
genome and sequencing. ``pairtools parse`` ignores shorter *gaps* and keeps
longer ones as "null" alignments. The maximal size of ignored *gaps* is set by
the ``--max-inter-align-gap`` flag (by default, 20bp).


Parse2
-------------------------

If your Hi-C has long reads, you may want to report all the alignments in the reads with ``pairtools parse2``.

The complex walks are DNA molecules containing more than one ligation junction that may end up in more than one alignment
on forward, reverse, or both reads:

.. figure:: _static/rescue_modes.svg
   :width: 60 %
   :alt: Different modes of reporting complex walks
   :align: center

   Different modes of reporting complex walks

``pairtools parse2`` detects such molecules and **rescues** them.

Briefly, the algorithm of complex ligation walks rescue detects all the unique ligation junctions, and do not report
the same junction as a pair multiple times. Importantly, these duplicated pairs might arise when both forward and reverse
reads read through the same ligation junction. However, these cases are successfully merged by ``pairtools parse2``:

.. figure:: _static/rescue_modes_readthrough.svg
   :width: 60 %
   :alt: Reporting complex walks in case of readthrough
   :align: center

   Reporing complex walks in case of readthrough

To restore the sequence of ligation events, there is a special field ``junction_index`` that can be reported as
a separate column of .pair file by setting ``--add-junction-index``. This field contains information on:

- the order of the junction in the recovered walk, starting from 5'-end of forward read
- type of the junction:

  - "u" - unconfirmed junction, right and left alignments in the pair originate from different reads (forward or reverse). This might be indirect ligation (mediated by other DNA fragments).
  - "f" - pair originates from the forward read. This is direct ligation.
  - "r" - pair originated from the reverse read. Direct ligation.
  - "b" - pair was sequenced at both forward and reverse read. Direct ligation.
With this information, the whole sequence of ligation events can be restored from the .pair file.


.. _section-single-ligation-rescue:

.. [1] Following the lead of `C-walks <https://www.nature.com/articles/nature20158>`_

.. [2] This procedure was first introduced in `HiC-Pro <https://github.com/nservant/HiC-Pro>`_ 
   and the in `Juicer <https://github.com/theaidenlab/juicer>`_ .
