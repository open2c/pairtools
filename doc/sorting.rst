Sorting pairs
=============

In order to enable efficient random access to Hi-C pairs, we **flip** and **sort** pairs. 
After sorting, interactions become arranged in the order of their genomic position, 
such that, for any given pair of regions, we easily find and extract all of their interactions.
And, after flipping, all artificially duplicated molecules (either during PCR or
in optical sequencing) end up in adjacent rows in sorted lists of interactions,
such that we can easily identify and remove them.

Sorting
-------

``pairtools sort`` arrange pairs in the order of (chrom1, chrom2, pos1, pos2).
This order is also known as *block sorting*, because all pairs between
any given pair of chromosomes become grouped into one continuous block.
Additionally, ``pairtools sort`` also sorts pairs with identical positions by
`pair_type`. This does not really do much for mapped reads, but it nicely splits
unmapped reads into blocks of null-mapped and multi-mapped reads.

We note that there is an alternative to block sorting, called *row sorting*, 
where pairs are sorted by (chrom1, pos1, chrom2, pos2). 
In `pairtools sort`, we prefer block-sorting since it cleanly separates cis 
interactions from trans ones and thus is a more optimal solution for typical
use cases.


Flipping
--------

In a typical paired-end experiment, *side1* and *side2* of a DNA molecule are
defined by the order in which they got sequenced.
Since this order is essentially random, any given Hi-C pair, e.g. 
(chr1, 1.1Mb; chr2, 2.1Mb), may appear in a reversed orientation, i.e.
(chr2, 2.1Mb; chr1, 1.1Mb). If we were to preserve this order of sides, interactions
between same loci would appear in two different locations of the sorted pair list,
which would complicate finding PCR/optical duplicates.

To ensure that Hi-C pairs with similar coordinates end up in the same location of the sorted list,
we **flip** pairs, i.e. we choose *side1* as the side with the lowest genomic coordinate. 
Thus, after flipping, for *trans* pairs (chrom1!=chrom2), order(chrom1)<order(chrom2);
and for *cis* pairs (chrom1==chrom2), pos1<=pos2.
In a matrix representation, flipping is equal to reflecting the lower triangle
of the Hi-C matrix onto its upper triangle, such that the resulting matrix 
is upper-triangular.

In `pairtools`, flipping is done during parsing - that's why ``pairtools parse``
requires a .chromsizes file that specifies the order of chromosomes for flipping.
Importantly, ``pairtools parse`` also flips one-sided pairs such that
side1 is always unmapped; and unmapped pairs such that side1 always has a "poorer"
mapping type (i.e. null-mapping<multi-mapping).


Chromosomal order for sorting and flipping
------------------------------------------

Importantly, the order of chromosomes for sorting and flipping can be different.
Specifically, ``pairtools sort`` uses the lexicographic order for chromosomes
(chr1, chr10, chr11, ..., chr2, chr21,...) instead of the "natural" order
(chr1, chr2, chr3, ...); at the same time, flipping is done in
``pairsamtools parse`` using the chromosomal order specified by the user.

``pairtools sort`` uses the lexicographic order for sorting chromosomes.
This order is used universally to sorting strings in all languages and tools [1]_, 
which makes it easy to design tools for indexing and searching in sorted pair lists.

At the same time, ``pairtools parse`` uses a custom user-provided chromosomal
order to flip pairs. For performance considerations, for flipping, we recommend
ordering chromosomes in a way that will be used in plotting contact maps.

.. [1] Unfortunately, many existing genomes use rather unconventional choices
   in chromosomal naming schemes. For example, in sacCer3, chromosomes are
   enumerated with Roman numerals; in dm3, big autosomes are split 4 different
   contigs each. Thus, it is impossible to design a universal algorithm that
   would order chromosomes in a "meaningful" way for all existing genomes.


