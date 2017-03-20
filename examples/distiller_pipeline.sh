#!/usr/bin/env bash
if [ $# -le 3 ] ; then
    echo "Usage: bash distiller_pipeline BWA_INDEX FASTQ_1 FASTQ_2 OUTPUT_PREFIX"
    echo ""
    echo "A example of a bash pipeline to align the sequencing data from a "
    echo "single Hi-C experiment."
    echo ""
    echo "positional arguments:"
    echo ""
    echo "BWA_INDEX       The path to a bwa index of the reference genome."
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
FASTQ1=$2
FASTQ2=$3
OUTPREFIX=$4

N_THREADS=8

UNMAPPED_SAM_PATH=${OUTPREFIX}.unmapped.bam
UNMAPPED_PAIRS_PATH=${OUTPREFIX}.unmapped.pairs.gz
NODUPS_SAM_PATH=${OUTPREFIX}.nodups.bam
NODUPS_PAIRS_PATH=${OUTPREFIX}.nodups.pairs.gz
DUPS_SAM_PATH=${OUTPREFIX}.dups.bam
DUPS_PAIRS_PATH=${OUTPREFIX}.dups.pairs.gz

DISTILLER_DIR="$(dirname ${BASH_SOURCE[0]})"
UTILS_DIR=${DISTILLER_DIR}/utils

bwa mem -SP -t "${N_THREADS}" "${INDEX}" "${FASTQ1}" "${FASTQ2}" | {
    # Classify Hi-C molecules as unmapped/single-sided/multimapped/chimeric/etc
    # and output one line per read, containing the following, separated by \\v:
    #  * triu-flipped pairs
    #  * read id
    #  * type of a Hi-C molecule
    #  * corresponding sam entries
    python ${UTILS_DIR}/sam_to_pairsam.py 
} | {
    # Block-sort pairs together with SAM entries
    bash ${UTILS_DIR}/pairsam_sort.sh
} | {
    # Set unmapped and ambiguous reads aside
    python ${UTILS_DIR}/pairsam_select.py pair_type CX,LL \
        --output-rest >( python ${UTILS_DIR}/pairsam_split.py \
            --output-pairs ${UNMAPPED_PAIRS_PATH} \
            --output-sam ${UNMAPPED_SAM_PATH} ) 
} | {
    # Remove duplicates
    python ${UTILS_DIR}/pairs_dedup.py \
        --output \
            >( python ${UTILS_DIR}/pairsam_split.py \
                --output-pairs ${NODUPS_PAIRS_PATH} \
                --output-sam ${NODUPS_SAM_PATH} ) \
        --output-dups \
            >( python ${UTILS_DIR}/pairsam_markasdup.py \
                | python ${UTILS_DIR}/pairsam_split.py \
                    --output-pairs ${DUPS_PAIRS_PATH} \
                    --output-sam ${DUPS_SAM_PATH} )
}

