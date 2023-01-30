`Pairtools stats` as a source of quality control metrics
===========================================

Overview
--------




Estimating library complexity
----------------------------

Pairtools assumes that each sequencing read is randomly chosen with
replacement from a finite pool of fragments in DNA library [1]_ [2]_.
With each new sequenced molecule, the expected number of observed unique molecules
increases according to a simple equation:

$$ U(N+1) = U(N) + (1-\frac{U(N)}{C}), $$

where $N$ is the number of sequenced molecules, $U(N)$ is the expected number
of observed unique molecules after sequencing $N$ molecules, and C is the library complexity.
This differential equation yields [1, 2]:

$$ \frac{U(N)}{C} = 1 - exp(-\frac{N}{C}), $$

which can be solved as

$$ C = \Re(lambert W(-\frac{exp(-\frac{1}{u})}{u}) + \frac{1}{u} $$

Library complexity can guide in the choice of sequencing depth of the library
and provide an estimate of library quality.


Optical duplicates
-----------------

Importantly, researchers can estimate the complexity of Hi-C libraries using only small QC
samples to decide if their quality permits deeper sequencing [library complexity papers].
These estimates, however, can be significantly biased by the presence of “optical” or
“clustering” duplicates. Such duplicates occur when a DNA sequencer erroneously splits
a signal from a single sequenced molecule into two; alternatively, a molecule located between
two adjacent tiles of a flowcell can be imaged twice, in both of the tiles [ref].

The rate of optical duplication depends on the technology and the operating conditions,
but not on the library complexity and sequencing depth. Thus, in small sequencing samples
optical duplication can severely inflate the observed levels of duplication,
resulting in underestimation of the library complexity.

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


.. [1] Picard. http://broadinstitute.github.io/picard/

.. [2] Thread: [Samtools-help] Pickard estimate for the size of a library - wrong or non-transparent? https://sourceforge.net/p/samtools/mailman/samtools-help/thread/DUB405-EAS154589A1ACEF2BE4C573D4592180@phx.gbl/
