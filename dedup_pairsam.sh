[ $# -ge 1 -a -f "$1" ] && INPUT="$1" || INPUT="-"

UNMAPPED_PATH=./test.unmapped.pairsam
DEDUP_PATH=./test.dedup.pairsam
DUPS_PATH=./test.dups.pairsam

AWK_SPLIT_UNMAPPED='
    BEGIN { FS="\v" }    
    {
        if ( $0 ~ /^#/ ) {print $0; print > unmappedpath} 
        else if ( $8 ~ /^CX|^LL/ ) {print } 
        else {print > unmappedpath}
    }
'

cat $INPUT | awk -v unmappedpath="$UNMAPPED_PATH" "$AWK_SPLIT_UNMAPPED" | \
    python dedup.py  \
    --out "$DEDUP_PATH" \
    --dupfile >(python tag_duplicates.py  > "$DUPS_PATH")

