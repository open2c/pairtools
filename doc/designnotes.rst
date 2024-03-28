Design notes
=============

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
- `pairtools` rely on a text format, .pairs, instead of hdf5/parquet-based
  tables or custom binaries. We went with a text format for a few reasons:

  - text tables enable easy access to data from any language and any tool. 
    This is especially important at the level of Hi-C pairs, the "rawest"
    format of information from a Hi-C experiment.
  - hdf5 and parquet have a few shortcomings that hinder their immediate use 
    in `pairtools`. Specifically, hdf5 cannot compress variable-length strings
    (which are, in turn, required to store sam alignments and some optional
    information on pairs) and parquet cannot append columns to existing files,
    modify datasets in place or store multiple tables in one file (which is
    required to keep table indices in the same file with pairs).
  - text tables have a set of well-developed and highly-optimized tools for
    sorting (Unix sort), compression (bgzip/lz4) and random access (tabix).
  - text formats enable easy streaming between individual command-line tools.
  
  Having said that, text formats have many downsides - they are bulky when
  not compressed, compression and parsing requires extra computational 
  resources, they cannot be modified in place and random access requires extra
  tools. In the future, we plan to develop a binary format based on existing
  container formats, which would mitigate these downsides.


CLI
---

- many `pairtools` perform multiple actions at once, which contradicts the
  "do one thing" philosophy of Unix command line. We packed multiple (albeit,
  related) functions into one tool to improve the performance of `pairtools`.
  Specifically, given the large size of Hi-C data, a significant fraction of time
  is spent on compression/decompression, parsing, loading data into memory and 
  sending it over network (for cloud/clusters). Packing multiple functions
  into one tool cuts down the amount of such time consuming operations.
- ``pairtools parse`` requires a .chromsizes file to know the order of chromosomes
  and perform pair flipping.
- `pairtools` use `bgzip <http://www.htslib.org/doc/bgzip.html>`_ compression by
  default instead of gzip. Using `bgzip` allows us to create an index with 
  `pairix <https://github.com/4dn-dcic/pairix>`_ and get random access to data.
- `paritools` have an option to compress outputs with
  `lz4 <https://en.wikipedia.org/wiki/LZ4_(compression_algorithm)>`_.
  `Lz4 is much faster and only slighly less efficient than gzip
  <https://catchchallenger.first-world.info/wiki/Quick_Benchmark:_Gzip_vs_Bzip2_vs_LZMA_vs_XZ_vs_LZ4_vs_LZO>`_.
  This makes lz4 a better choice for passing data between individual pairtools
  before producing final result (which, in turn, requires bgzip compression).


