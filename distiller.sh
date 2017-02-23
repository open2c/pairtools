INDEX=$1                                                   
FASTQ1=$2                                                  
FASTQ2=$3                                                  
OUTNAME=$4                                                 

~/programming/bwa/bwa mem -SP $INDEX $FASTQ1 $FASTQ2 | {
# drop the supplementary alignments b/c the information is already contained in the primary one
     samtools view -h -F 2048
} | {
# save unmapped/single-sided/multimapped/abnormal chimeras 
# and output flipped pairs followed by the two bam entries, separated by \v
    python split_sams_append_pairs.py 
        --header $OUTNAME.header.sam
        --unmapped >(samtools view -bS - > $OUTNAME.unmapped.bam)
        --singlesided >(samtools view -bS - > $OUTNAME.singlesided.bam)
        --multimapped >(samtools view -bS - > $OUTNAME.mutimapped.bam)
        --abnormal-chimera >(samtools view -bS - >$OUTNAME.abnormal_chimera.bam)
        
} | {
# sort pairs together with bams
    sort -k 1,1 -k 4,4 -k 2,2n -k 5,5n --field-separator='\v' 
} | {
# remove duplicates 
    python dedup.py 
        --out >(python split_pairs.py 
                --header $OUTNAME.header.sam
                --out-pairs $OUTNAME.dedup.pairs
                --out-sam >(samtools view -bS - > $OUTNAME.dedup.bam)
        )
        --dupfile >(python split_pairs.py
                --header $OUTNAME.header.sam
                --out-pairs /dev/null
                --out-sam >(samtools view -bS - > $OUTNAME.dups.bam)
        )
}

