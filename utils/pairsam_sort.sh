#!/usr/bin/env bash
if [ $# -ge 1 -a -f "$1" ] ; then
    DIR="$(dirname "${BASH_SOURCE[0]}")"
    { bash $DIR/pairsam_get_header.sh $1 ; bash $DIR/pairsam_skip_header.sh $1 | sort -k 1,1 -k 4,4 -k 2,2n -k 5,5n -k 8,8 --field-separator=$'\v' ; }
else 
    awk '/^#/ {print $0;next} {print $0 | "sort -k 1,1 -k 4,4 -k 2,2n -k 5,5n -k 8,8 --field-separator=\v"}'
fi 

