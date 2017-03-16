[ $# -ge 1 -a -f "$1" ] && INPUT="$1" || INPUT="-"

UNMAPPED_PATH=./test.unmapped.pairsam
DEDUP_PATH=./test.dedup.pairsam
DUPS_PATH=./test.dups.pairsam

AWK_SPLIT_UNMAPPED='{
    if ( $0 ~ /^#/ ) {print $0; print > unmappedpath} 
    else if ( $0 ~ /^CX|^LL/ ) {print } 
    else {print > unmappedpath}
}'

cat $INPUT | awk -v unmappedpath="$UNMAPPED_PATH" "$AWK_SPLIT_UNMAPPED" | \
    python dedup.py --c1 1 --c2 4 --p1 2 --p2 5 --s1 3 --s2 6  \
    --out "$DEDUP_PATH" \
    --dupfile >(python tag_duplicates.py  > "$DUPS_PATH")

