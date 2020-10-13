Quickstart
==========

Install `pairtools` and all of its dependencies using the 
`conda <https://conda.io/docs/user-guide/install/index.html>`_ package manager and 
the `bioconda <https://bioconda.github.io/index.html>`_ channel for bioinformatics 
software.

.. code-block:: bash

    $ conda install -c conda-forge -c bioconda pairtools

Setup a new test folder and download a small Hi-C dataset mapped to sacCer3 genome:

.. code-block:: bash

    $ mkdir /tmp/test-pairtools
    $ cd /tmp/test-pairtools
    $ wget https://github.com/open2c/distiller-test-data/raw/master/bam/MATalpha_R1.bam

Additionally, we will need a .chromsizes file, a TAB-separated plain text table describing the names, sizes and the order of chromosomes in the genome assembly used during mapping:

.. code-block:: bash

    $ wget https://raw.githubusercontent.com/open2c/distiller-test-data/master/genome/sacCer3.reduced.chrom.sizes

With `pairtools parse`, we can convert paired-end sequence alignments stored in .sam/.bam format into .pairs, a TAB-separated table of Hi-C ligation junctions:

.. code-block:: bash

    $ pairtools parse -c sacCer3.reduced.chrom.sizes -o MATalpha_R1.pairs.gz --drop-sam MATalpha_R1.bam 

Inspect the resulting table:

.. code-block:: bash

    $ less MATalpha_R1.pairs.gz

