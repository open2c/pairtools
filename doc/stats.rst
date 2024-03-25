`Pairtools stats` as a source of quality control metrics
===========================================

Overview
--------

`pairtools stats` produces a human-readable nested dictionary of statistics stored in
a YAML file or a tab-separated text table (specified through the parameters).

When calculating statistics, any number of filters can be applied to generate separate
statistics for different categories of pairs, for example they can be filtered by the
read mapping quality (mapq values). These are then stored as separate sections of the
output file.

- **Global statistics** include:
    - number of pairs (total, unmapped, single-side mapped, etc.),
    - total number of different pair types (UU, NN, NU, and others, see `Pair types in pairtools docs <https://pairtools.readthedocs.io/en/latest/formats.html#pair-types>`_),
    - number of contacts between all chromosome pairs

- **Summary statistics** include:
    - fraction of duplicates
    - fraction of cis interactions (at different minimal distance cutoffs) out of total
    - estimation of library complexity

Summary statistics can inform you about the quality of the data.
For example, more trans interactions can be a sign of problems with the 3C+ procedure and lower signal-to-noise ratio.
Substantial mapping to mitochondrial chromosome (chrM) might be a sign of random ligation.

- **P(s), or scaling.**  The dependence of contact frequency on the genomic
distance referred to as the P(s) curve or scaling, which is a rich source of both biologically relevant information and technical quality of 3C+ experiments.
The shape of P(s) is often used to characterize mechanisms of genome folding and reveal issues with QC.

Interactive visualization of stats with MultiQC
---------

Install `multiqc`:

.. code-block:: bash

    pip install --upgrade --force-reinstall git+https://github.com/open2c/MultiQC.git

Note that (for now) the pairtools module for MultiQC is only available in the open2C fork and not in the main MultiQC repository.

Run MultiQC in a folder with one or multiple .stats files:

.. code-block:: bash

    multiqc .


This will produce a nice .html file with interactive graphical summaries of the stats.


Estimating library complexity
----------------------------

Pairtools assumes that each sequencing read is randomly chosen with
replacement from a finite pool of fragments in DNA library [1]_ [2]_.
With each new sequenced molecule, the expected number of observed unique molecules
increases according to a simple equation:

.. math::

    U(N+1) = U(N) + \left(1 - \frac{U(N)}{C} \right),

where :math:`N` is the number of sequenced molecules, :math:`U(N)` is the expected number
of observed unique molecules after sequencing :math:`N` molecules, and :math:`C` is the library complexity.
This differential equation yields [1, 2]:

.. math::
    
    {U(N) \over C} = 1 - exp\left( - \frac{N}{C} \right),

which can be solved as

.. math::

    C = \Re \left( W_{Lambert} \left( - \frac{ \exp\left( - \frac{1}{U} \right) } {U} \right) \right) + \frac{1}{U}

Library complexity can guide in the choice of sequencing depth of the library
and provide an estimate of library quality.


Illumina sequencing duplicates
-----------------

Importantly, you can estimate the complexity of Hi-C libraries using only small QC
samples to decide if their quality permits deeper sequencing [3]_.
These estimates, however, can be significantly biased by the presence of “optical” or
“clustering” duplicates. Such duplicates occur as artefacts of the sequencing procedure.
Optical duplicates appear in data generated on sequencers with non-patterned flowcells in
cases the instrument either erroneously splits a signal from a single sequenced molecule
into two. On the other hand, clustering duplicates appear on patterned flowcells, when
during cluster generation a cluster occupies adjacent nanowells. [4]_.

The rate of optical and clustering duplication depends on the technology and the operating
conditions (e.g. molarity of the library loaded onto the flowcell), but not on the
library complexity or sequencing depth. Thus, in small sequencing samples in particular
the clustering duplication on recent Illumina instruments can severely inflate the
observed levels of duplication [5]_, resulting in underestimation of the library complexity.

While the frequency of PCR duplicates increases with sequencing depth,
optical or clustering duplication levels may stay constant for a particular sequencer,
provided the library is loaded at the same molarity. This means that the high frequency of
clustering duplicates on the NovaSeq leads to severe underestimation of library complexity
in the pilot runs. In particular, the recent models of Illumina sequencers with patterned
flowcells (such as NovaSeq) suffer from increased clustering duplication rate, which may
far exceed the level of PCR duplication.

Luckily, optical and clustering duplicates can be distinguished from the PCR ones,
as the former are located next to each other on the sequencing flow cell.
In case of Illumina sequencers, pairtools dedup can infer the positions of sequencing
reads from their IDs and focuses on geometrically distant duplicates to produce unbiased
estimates of PCR duplication and library complexity.  Although SRA does not store original
read IDs from the sequencer, this analysis is possible when pairtools is run on a dataset
with original Illumina-generated read IDs.
Note that in our experience even when accounting for optical/clustering duplicates, the
complexity can be greatly underestimated, but is still a useful measurement to choose the
most complex libraries.


.. [1] Picard. http://broadinstitute.github.io/picard/

.. [2] Thread: [Samtools-help] Pickard estimate for the size of a library - wrong or non-transparent? https://sourceforge.net/p/samtools/mailman/samtools-help/thread/DUB405-EAS154589A1ACEF2BE4C573D4592180@phx.gbl/

.. [3] Rao, S. S. P. et al. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. Cell 159, 1665–1680 (2014).

.. [4] Duplicates on Illumina. BioStars. https://www.biostars.org/p/229842/
.. [5] Illumina Patterned Flow Cells Generate Duplicated Sequences. https://sequencing.qcfail.com/articles/illumina-patterned-flow-cells-generate-duplicated-sequences/