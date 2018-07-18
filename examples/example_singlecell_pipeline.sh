#!/usr/bin/env bash
if [ $# -le 3 ] ; then
    echo "Usage: bash example_pipeline.sh BWA_INDEX FASTQ_1 FASTQ_2 OUTPUT_PREFIX"
    echo ""
    echo "A example of a bash pipeline to align the sequencing data from a "
    echo "single Hi-C experiment."
    echo ""
    echo "positional arguments:"
    echo ""
    echo "BWA_INDEX       The path to a bwa index of the reference genome."
    echo "CHROM_SIZES     The path to a file with chromosome sizes."
    echo "FASTQ_1         The path to a fastq file with the sequences of "
    echo "                the first side of Hi-C molecules."
    echo "FASTQ_2         The path to a fastq file with the sequences of "
    echo "                the second side of Hi-C molecules."
    echo "OUTPUT_PREFIX   The prefix to the paths of generated outputs. "
    echo ""
    echo ""

    exit 0
fi

set -o errexit
set -o nounset
set -o pipefail

INDEX=$1
CHROM_SIZES=$2
FASTQ1=$3
FASTQ2=$4
OUTPREFIX=$5

N_THREADS=8

UNMAPPED_SAM_PATH=${OUTPREFIX}.unmapped.bam
UNMAPPED_PAIRS_PATH=${OUTPREFIX}.unmapped.pairs.gz
NODUPS_SAM_PATH=${OUTPREFIX}.nodups.bam
NODUPS_PAIRS_PATH=${OUTPREFIX}.nodups.pairs.gz
DUPS_SAM_PATH=${OUTPREFIX}.dups.bam
DUPS_PAIRS_PATH=${OUTPREFIX}.dups.pairs.gz
LOWFREQPAIRS_SAM_PATH=${OUTPREFIX}.lowfreq.bam
LOWFREQPAIRS_PAIRS_PATH=${OUTPREFIX}.lowfreq.pairs.gz
HIGHFREQPAIRS_SAM_PATH=${OUTPREFIX}.highfreq.bam
HIGHFREQPAIRS_PAIRS_PATH=${OUTPREFIX}.highfreq.pairs.gz

bwa mem -SP -t "${N_THREADS}" "${INDEX}" "${FASTQ1}" "${FASTQ2}" | {
    # Classify Hi-C molecules as unmapped/single-sided/multimapped/chimeric/etc
    # and output one line per read, containing the following, separated by \\v:
    #  * triu-flipped pairs
    #  * read id
    #  * type of a Hi-C molecule
    #  * corresponding sam entries
    pairtools parse "{CHROM_SIZES}"
} | {
    # Block-sort pairs together with SAM entries
    pairtools sort
} | {
    # Set unmapped and ambiguous reads aside
    pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
        --output-rest >( pairtools split \
            --output-pairs ${UNMAPPED_PAIRS_PATH} \
            --output-sam ${UNMAPPED_SAM_PATH} ) 
} | {
    # Remove duplicates
    pairtools dedup \
        --output-dups \
            >( pairtools markasdup \
                | pairtools split \
                    --output-pairs ${DUPS_PAIRS_PATH} \
                    --output-sam ${DUPS_SAM_PATH} )
} | {
    # Remove high frequency interactors
    pairtools multifilter \ 
        --output \
            >( pairtools split \
                --output-pairs ${LOWFREQ_PAIRS_PATH} \
                --output-sam ${LOWFREQ_SAM_PATH} ) \
        --output-high-frequency-interactors \
            >( pairtools markasdup \
                | pairtools split \
                    --output-pairs ${HIGHFREQPAIRS_PAIRS_PATH} \
                    --output-sam ${HIGHFREQPAIRS_SAM_PATH} )
}

