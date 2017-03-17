#!/bin/bash
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


bwa mem -SP "$INDEX" "$FASTQ1" "$FASTQ2" | {
# classify Hi-C molecules as unmapped/single-sided/multimapped/chimeric/etc
# and output one line per read, containing flipped pairs, type of a Hi-C molecule
# and the corresponding sam entries, all separated by \v
    python classify_reads.py 
} | {
# lexicographic block-sort pairs together with sam entries
    bash pairsam_sort.sh
} | {
# remove duplicates and split pairs and sams
    cat - \
        | bash pairsam_divertunmapped.sh \
            >( python pairsam_split.py \
                --out-pairs ${UNMAPPED_PAIRS_PATH} \
                --out-sam ${UNMAPPED_SAM_PATH} ) \
        | python dedup_pairs.py  \
            --out \
                >( python pairsam_split.py \
                    --out-pairs ${NODUPS_PAIRS_PATH} \
                    --out-sam ${NODUPS_SAM_PATH} )\
            --dupfile \
                >( python pairsam_markasdup.py | \
                    python pairsam_split.py \
                     --out-pairs ${DUPS_PAIRS_PATH} \
                     --out-sam ${DUPS_SAM_PATH} )
}
