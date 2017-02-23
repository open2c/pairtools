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
    python filter_addpairs.py 
        --header $OUTNAME.header.sam
        --unmapped $OUTNAME.unmapped.bam
        --singlesided $OUTNAME.singlesided.bam
        --multimapped $OUTNAME.mutimapped.bam
        --abnormal-chimera $OUTNAME.abnormal_chimera.bam
        
} #| {
# sort pairs together with bams
#    | sort -k 1,1 -k 2,2n -k 4,4 -k 5,5n --field-separator='\v' \

# remove duplicates in the pair list and flag duplicates in bams
#    | dedup  \
# split off pairs into a separate list
#    | python tee_contacts.py ${OUTNAME}.validPairs.txt.gz \
#    | samtools -bS - > ${OUTNAME}.validPairs.bam

