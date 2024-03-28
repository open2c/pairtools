Workflows and Parameters
========================

This page provides guidance on using pairtools for the most common Hi-C protocols and 
helps users fine-tune the pipeline for different variations of the Hi-C protocol. 
It covers recommended parameters and best practices for processing Hi-C data using pairtools.

Typical Hi-C Workflow
----------------------

A typical pairtools workflow for processing standard Hi-C data is outlined below. 
Please, note that this is a shorter version. For a detailed reproducible example, please, check the Jupyter notebook "Pairtools Walkthrough".

1. Align sequences to the reference genome with ``bwa mem``:
   
    .. code-block:: console

        bwa mem -SP index_file input.R1.fastq input.R2.fastq > input.sam

2. Parse alignments into Hi-C pairs using ``pairtools parse``:

    .. code-block:: console 

        pairtools parse -c /path/to/chrom_sizes -o output.pairs.gz input.sam

3. Sort pairs using ``pairtools sort``:


    .. code-block:: console


        pairtools sort --nproc 8 -o output.sorted.pairs.gz output.pairs.gz

4. Detect and remove duplicates using ``pairtools dedup`` and generate statistics:

    .. code-block:: console

        pairtools dedup \
        --output output.nodups.pairs.gz \
        --output-dups output.dups.pairs.gz \
        --output-unmapped output.unmapped.pairs.gz 
        --output-stats output.dedup.stats \
        output.sorted.pairs.gz

5. Aggregate into a cooler file:

    .. code-block:: console

        cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 /path/to/chrom_sizes:1000 output.nodups.pairs.gz output.1000.cool




Recommended pairtools parameters for standard Hi-C protocols
------------------------------------------------------------

To adapt the standard workflow for common variations of the Hi-C protocol, consider adjusting the following parameters:

1. ``pairtools parse --walks-policy``: 

   This parameter determines how pairtools parse handles reads with multiple alignments (walks). We recommend specifying the value explicitly, as the default has changed between versions of ``pairtools parse``.
    
   Our current recommendation is to use ``--walks-policy 5unique``, which is the default setting in the latest version of pairtools. With this option, pairtools parse reports the two 5'-most unique alignments on each side of a paired read as a pair. 

   This option increases the number of reported pairs compared to the most conservative ``--walks-policy mask``. However, it's important to note that ``5unique`` can potentially report pairs of non-directly ligated fragments (i.e., two fragments separated by one or more other DNA fragments). Such non-direct (also known as "higher-order" or "nonadjacent") ligations have slightly different statistical properties than direct ligations, as illustrated in several Pore-C papers  [`1 <https://www.biorxiv.org/content/10.1101/833590v1.full>`_ , `2 <https://www.nature.com/articles/s41467-023-36899-x>`_].

   An alternative is the ``--walks-policy 3unique`` policy, which reports the two 3'-most unique alignments on each side of 
   a paired read as a pair, thus decreasing the chance of reporting non-direct ligations. 
   However, ``3unique`` may not work well in situations where the combined length of a read pair is longer than the length of a DNA fragment (e.g. long read experiments). 
   In this case, the 3' sides of the two reads will cover the same locations in the DNA molecule, and the 3' alignments may end up identical.
    
   Finally, the experimental ``--walks-policy all`` option reports all alignments of a read pair as separate pairs. 
   This option maximizes the number of reported pairs. 
   The downside is that it breaks the assumption that there is only one pair per read, 
   which is not compatible with retrieval of .sam records from .pairsam output and may also complicate the interpretation of pair statistics.

2. ``pairtools select "(mapq1>=30) and (mapq2>=30)"``: 

   This filtering command selects only pairs with high-quality alignments, 
   where both reads in a pair have a mapping quality (MAPQ) score of 30 or higher. 
   Applying this filter helps remove false alignments between partially homologous sequences, which often cause artificial high-frequency interactions in Hi-C maps. 
   This step is essential for generating maps for high-quality dot calls.

   Note that we recommend storing the most comprehensive, unfiltered list of pairs and applying the filter on the fly prior to contact aggregation:

    .. code-block:: console

        pairtools select "(mapq1>=30) and (mapq2>=30)" output.nodups.pairs.gz | \
            cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 chromsizes.txt:1000 - output.mapq_30.1000.cool


Technical tips
--------------

- **Pipe between commands to save space and I/O throughput**

    Use Unix pipes to connect the output of one command directly to the input of the next command in the pipeline. 
    This eliminates the need to store intermediate files on disk, saving storage space and reducing I/O overhead.
    Specifically, mapping, parsing, sorting and deduplication can all be connected into a single pipeline:

    .. code-block:: console

        bwa mem -SP index input.R1.fastq input.R2.fastq | \
        pairtools parse -c chromsizes.txt | \
        pairtools sort | \
        pairtools dedup | \
            --output output.nodups.pairs.gz \
            --output-dups output.dups.pairs.gz \
            --output-unmapped output.unmapped.pairs.gz 
            --output-stats output.dedup.stats

- **Use recommended compression for efficient storage and processing.** .sam, .pairs and .pairsam files are text-based format that are rather inefficient and slow to process.  
  Pairtools recognize .bam, .gz and .lz4 file extensions and automatically compress and decompress files on the fly.
  Compression saves space, and reduces I/O overhead at a relatively minor CPU cost.

- **Parallelize tasks and manage resources effectively for faster execution.**
  Each pairtool has the CLI flags --nproc-in and --nproc-out to control the number of cores dedicated 
  to input decompression and output compression. Additionally, `pairtools sort` parallelizes sorting with `--nproc`.ß

Advanced Workflows
------------------

For more advanced workflows, please check the following projects:

- `Distiller-nf <https://github.com/open2c/distiller-nf>`_ is a feature-rich Open2C Hi-C processing pipeline for the Nextflow workflow manager.
- `Distiller-sm <https://github.com/open2c/distiller-sm>`_ is a similarly feature-rich and optimized pipeline implemented in Snakemake.
