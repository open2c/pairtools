#!/bin/bash
PREFIX=$1

DIR="$(dirname "${BASH_SOURCE[0]}")"

UNMAPPED_SAM_PATH=${PREFIX}.unmapped.bam
UNMAPPED_PAIRS_PATH=${PREFIX}.unmapped.pairs.gz
NODUPS_SAM_PATH=${PREFIX}.nodups.bam
NODUPS_PAIRS_PATH=${PREFIX}.nodups.pairs.gz
DUPS_SAM_PATH=${PREFIX}.dups.bam
DUPS_PAIRS_PATH=${PREFIX}.dups.pairs.gz

cat - | \
    bash $DIR/split_unmapped.sh \
        >( python split_pairsam.py \
            --out-pairs ${UNMAPPED_PAIRS_PATH} \
            --out-sam ${UNMAPPED_SAM_PATH} ) | \
    python dedup.py  \
        --out \
            >( python split_pairsam.py \
                --out-pairs ${NODUPS_PAIRS_PATH} \
                --out-sam ${NODUPS_SAM_PATH} )\
        --dupfile \
            >( python tag_duplicates.py | python split_pairsam.py \
                --out-pairs ${DUPS_PAIRS_PATH} \
                --out-sam ${DUPS_SAM_PATH} )
