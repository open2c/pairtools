INDEX=$1                                                   
FASTQ1=$2                                                  
FASTQ2=$3                                                  
OUTNAME=$4                                                 

~/programming/bwa/bwa mem -SP $INDEX $FASTQ1 $FASTQ2 | {
# save unmapped/single-sided/multimapped/abnormal chimeras 
# and output lines with flipped pairs and the two sam entries, separated by \v
    python classify_sam.py \
        --header $OUTNAME.header.sam \
        --unmapped >(samtools view -bS - > $OUTNAME.unmapped.bam) \
        --singlesided >(samtools view -bS - > $OUTNAME.singlesided.bam) \
        --multimapped >(samtools view -bS - > $OUTNAME.mutimapped.bam) \
        --abnormal-chimera >(samtools view -bS - >$OUTNAME.abnormal_chimera.bam) 
} | {
# lexicographic block-sort pairs together with sam entries
    sort -k 1,1 -k 4,4 -k 2,2n -k 5,5n --field-separator=\v
} | {
# remove duplicates and split pairs and sams
    python dedup_pairsam.py \
        --out >(python split_pairsam.py \
                --header $OUTNAME.header.sam \
                --out-pairs $OUTNAME.dedup.pairs \
                --out-sam >(samtools view -bS - > $OUTNAME.dedup.bam)
        ) \
        --dupfile >(python split_pairsam.py \
                --header $OUTNAME.header.sam \
                --out-pairs /dev/null \
                --out-sam >(samtools view -bS - > $OUTNAME.dups.bam)
        )
}

