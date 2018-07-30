Technical notes
===============

Designing scientific software and formats requires making a multitude of 
tantalizing technical decisions and compromises. Often, the reasons behind a 
certain decision are non-trivial and convoluted, involving many factors.
Here, we collect the notes and observations made during the desing stage of 
`pairtools` and provide a justification for most non-trivial decisions.
We hope that this document will elucidate the design of `pairtools` and
may prove useful to developers in their future projects.

.pairs format
-------------

The motivation behind some of the technical decisions in the pairtools' flavor
of .pairs/.pairsam:

- `pairtools` can store SAM entries together with the Hi-C pair information in 
  .pairsam files. Storing pairs and alignments in the same row enables easy 
  tagging and filtering of paired-end alignments based on their Hi-C 
  information.
- `pairtools` use the exclamation mark "!" instead of '.' as 'chrom' of 
  unmapped reads because it has the lowest lexicographic sorting order among all
  characters. The use of '0' and '-' in the 'pos' and 'strand' fields of unmapped
  reads allows us to keep the types of these fields as 'unsigned int' and
  enum{'+','-'}, respectively.
- "rescued" pairs have two types "UR" and "RU" instead of just "RU". We chose
  this design because rescued pairs are two-sided and thus are flipped based on 
  (chrom, pos), and not based on the side types. With two pair types "RU" and "UR", 
  `pairtools` can keep track of which side of the pair was rescued.
- in "rescued" pairs, the type "R" is assigned to the non-chimeric side.
  This may seem counter-intuitive at first, since it is the chimeric side that
  gets rescued, but this way `pairtools` can keep track of the type of the
  5' alignment on the chimeric side (the alignment on the non-chimeric side
  has to be unique for the pair to be rescued).
- `pairtools` rely on a text format, .pairs, instead of a hdf5/parquet-based
  table or a custom binary. We went with a text format for a few reasons:

  - text tables enable easy access to data from any language and any tool. 
    This is especially important at the level of Hi-C pairs, the "rawest"
    format of information from a Hi-C experiment.
  - hdf5 and parquet have a few shortcomings that hinder their immediate use 
    in `pairtools`. Specifically, hdf5 cannot compress variable-length strings
    (which are, in turn, required to store sam alignments and some optional
    information on pairs) and parquet cannot append columns to existing files 
    and store multiple tables in one file (which is required to keep table 
    indices in the same file with pairs).
  - text tables have a set of well-developed and highly-optimized tools for
    sorting (Unix sort), compression (bgzip/lz4) and random access (tabix).
  - text formats allow streaming between individual tools.


Sorting
-------

- `pairtools` sort pairs by chromosomal blocks 
  (the sorting order is chrom1, chrom2, pos1, pos2) instead of rows (i.e. 
  chrom1, pos1, chrom2, pos2) Block-sorting  slower query on 
- `pairtools` sort chromosomes in the lexicographic order 
  (chr1, chr10, chr11, ..., chr2, chr21,...) instead of the "natural" order
  (chr1, chr2, chr3, ...). This is done 


CLI
---

- many `pairtools` perform multiple actions at once, which contradicts the
  "do one thing" philosophy of Unix command line. We had to lump multip
- `pairtools parse` requires a .chromsizes file to know the order of chromosomes
  and perform pair flipping.

- `pairtools` use `bgzip <http://www.htslib.org/doc/bgzip.html>`_ compression by
  default instead of gzip. 
- additionally, we enable 
  `lz4 <https://en.wikipedia.org/wiki/LZ4_(compression_algorithm)>`_ compression
  https://catchchallenger.first-world.info/wiki/Quick_Benchmark:_Gzip_vs_Bzip2_vs_LZMA_vs_XZ_vs_LZ4_vs_LZO

