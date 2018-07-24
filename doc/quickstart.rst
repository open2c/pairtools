Quickstart
==========

Install `pairtools` and all of its dependencies using the 
[conda](https://conda.io/docs/user-guide/install/index.html) package manager and 
the [bioconda](https://bioconda.github.io/index.html) channel for bioinformatics 
software.
```sh
$ conda install -c conda-forge -c bioconda pairtools
```

Setup a new test folder and download a small Hi-C dataset mapped to sacCer3 genome:
```bash
$ mkdir /tmp/test-pairtools
$ cd /tmp/test-pairtools
$ wget https://github.com/mirnylab/distiller-test-data/raw/master/bam/MATalpha_R1.bam
```

Additionally, we will need a .chromsizes file, a TAB-separated plain text table describing the names, sizes and the order of chromosomes in the genome assembly used during mapping:
```bash
$ wget https://raw.githubusercontent.com/mirnylab/distiller-test-data/master/genome/sacCer3.reduced.chrom.sizes
```

With `pairtools parse`, we can convert paired-end sequence alignments stored in .sam/.bam format into .pairs, a TAB-separated table of Hi-C ligation junctions:

```bash
$ pairtools parse -c sacCer3.reduced.chrom.sizes -o MATalpha_R1.pairs.gz --drop-sam MATalpha_R1.bam 
```

Inspect the resulting table:

```bash
$ less MATalpha_R1.pairs.gz
```

