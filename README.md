# distiller

distiller is a simple command-line framework to process sequencing information 
from a Hi-C experiment.

distiller provides tools to perform the following operations:
- map the sequences of Hi-C molecules to reference genomes
- filter mapped sequences for valid Hi-C pairs
- remove PCR duplicates 
- tag .sam entries with the Hi-C specific information
- produce sorted lists of pairs in a form convenient for downstream analyses

distiller provides a simple bash pipeline.

