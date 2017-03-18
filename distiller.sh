#!/usr/bin/env bash
set -o errexit
set -o nounset
set -o pipefail
INDEX=$1
FASTQ1=$2
FASTQ2=$3
OUTPREFIX=$4
UNMAPPED_SAM_PATH=${OUTPREFIX}.unmapped.bam
UNMAPPED_PAIRS_PATH=${OUTPREFIX}.unmapped.pairs.gz
NODUPS_SAM_PATH=${OUTPREFIX}.nodups.bam
NODUPS_PAIRS_PATH=${OUTPREFIX}.nodups.pairs.gz
DUPS_SAM_PATH=${OUTPREFIX}.dups.bam
DUPS_PAIRS_PATH=${OUTPREFIX}.dups.pairs.gz


bwa mem -SP "${INDEX}" "${FASTQ1}" "${FASTQ2}" | {
    # Classify Hi-C molecules as unmapped/single-sided/multimapped/chimeric/etc
    # and output one line per read, containing the following, separated by \\v:
    #  * triu-flipped pairs
    #  * type of a Hi-C molecule
    #  * corresponding sam entries
    python sam_to_pairsam.py 
} | {
    # Block-sort pairs together with SAM entries
    bash pairsam_sort.sh
} | {
    # Set unmapped and ambiguous reads aside
    bash pairsam_divertunmapped.sh \
        >( python pairsam_split.py \
            --out-pairs ${UNMAPPED_PAIRS_PATH} \
            --out-sam ${UNMAPPED_SAM_PATH} ) \
} | {
    # Remove duplicates
    python pairs_dedup.py \
        --out \
            >( python pairsam_split.py \
                --out-pairs ${NODUPS_PAIRS_PATH} \
                --out-sam ${NODUPS_SAM_PATH} ) \
        --dupfile \
            >( python pairsam_markasdup.py \
                | python pairsam_split.py \
                    --out-pairs ${DUPS_PAIRS_PATH} \
                    --out-sam ${DUPS_SAM_PATH} )
}
