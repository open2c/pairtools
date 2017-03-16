INDEX=$1                                                   
FASTQ1=$2                                                  
FASTQ2=$3                                                  
OUTNAME=$4

~/programming/bwa/bwa mem -SP $INDEX $FASTQ1 $FASTQ2 | {
# save unmapped/single-sided/multimapped/abnormal chimeras 
# and output lines with flipped pairs and the two sam entries, separated by \v
    python classify_sam.py 
} | {
# lexicographic block-sort pairs together with sam entries
    bash sort_pairsam.sh
} | {
# remove duplicates and split pairs and sams
    bash dedup_split_pairsam.sh
}

