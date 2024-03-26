Workflows and Parameters
========================

Introduction
------------

This page provides guidance on using pairtools for the most common Hi-C protocols and 
helps users fine-tune the pipeline for different variations of the Hi-C protocol. 
It covers recommended parameters and best practices for processing Hi-C data using pairtools.

Standard Hi-C Workflow
----------------------

A typical pairtools workflow for processing standard Hi-C data is outlined below. 
Please, note that this is a shorter version; you can find a more detailed and reproducible example in chapter :ref:`./examples/pairtools_walkthrough.ipynb`.

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

Together, these steps can be stringed into a simple two-step pipeline:

   .. code-block:: console
        bwa mem -SP index input.R1.fastq input.R2.fastq | \
        pairtools parse -c chromsizes.txt | \
        pairtools sort | \
            --output output.nodups.pairs.gz \
            --output-dups output.dups.pairs.gz \
            --output-unmapped output.unmapped.pairs.gz 
            --output-stats output.dedup.stats
        cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 chromsizes.txt:1000 output.nodups.pairs.gz output.1000.cool


Variations of Standard Hi-C
---------------------------

To adapt the standard workflow for common variations of the Hi-C protocol, consider adjusting the following parameters:

1. ``--min-mapq`` and ``--max-mapq``: These parameters control the minimum and maximum mapping quality (MAPQ) thresholds for accepting alignments. For example:

   - Increase ``--min-mapq`` to reduce noise from low-quality alignments in protocols with lower library complexity.
   - Decrease ``--max-mapq`` to retain multi-mapping reads in protocols where such reads are informative, such as in repetitive regions.

2. ``--walks-policy``: This parameter determines how ``pairtools parse`` handles reads with multiple alignments (walks). For example:

   - Use ``--walks-policy all`` to keep all walks for protocols with lower library complexity or when analyzing repetitive regions.

3. ``--max-inter-align-gap``: This parameter sets the maximum acceptable gap between two alignments in a walk. For example:

   - Increase ``--max-inter-align-gap`` to retain more informative walks in protocols with longer fragments or reads.

4. ``--no-rescue-walks``: This flag disables the rescue of truncated walks. For example:

   - Use ``--no-rescue-walks`` to reduce false-positive pairs in protocols with high levels of chimeric reads or ligation errors.

Other Hi-C Protocols
--------------------

Pairtools can be used for non-standard Hi-C protocols with additional tools and parameters:

- Single-cell Hi-C: Use ``pairtools phase`` for data mapped to diploid genomes.
- Micro-C: Adjust ``pairtools dedup`` settings for different amplification strategies.
- HiChIP: Incorporate ``pairtools restrict`` for enzyme-specific analysis.

Best Practices and Tips
-----------------------

- Use recommended file formats and compression for efficient storage and processing.
- Parallelize tasks and manage resources effectively for faster execution.
- Troubleshoot common issues by referring to the documentation and seeking help from the community.

Example Workflows
-----------------

Example workflows for common Hi-C data processing scenarios are available in the `examples` directory of the pairtools repository. Each example includes sample datasets, step-by-step instructions, and example output files.

FAQs
----

For frequently asked questions related to using pairtools for different Hi-C protocols and experimental designs, please refer to the FAQ section of the documentation.