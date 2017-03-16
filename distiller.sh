INDEX=$1                                                   
FASTQ1=$2                                                  
FASTQ2=$3                                                  
OUTPREFIX=$4

~/programming/bwa/bwa mem -SP "$INDEX" "$FASTQ1" "$FASTQ2" | {
# classify Hi-C molecules as unmapped/single-sided/multimapped/chimeric/etc
# and output one line per read, containing flipped pairs, type of a Hi-C molecule
# and the corresponding sam entries, all separated by \v
    python classify_sam.py 
} | {
# lexicographic block-sort pairs together with sam entries
    bash sort_pairsam.sh
} | {
# remove duplicates and split pairs and sams
    bash dedup_split_pairsam.sh "$OUTPREFIX"
}

