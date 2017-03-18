#!/usr/bin/env bash
AWK_SPLIT_UNMAPPED_SCRIPT='
    BEGIN { FS="\v" } 
    {
        if ( $0 ~ /^#/ ) {print $0; print > unmappedpath} 
        else if ( $8 ~ /CX|LL/ ) {print } 
        else {print > unmappedpath}
    }
'

awk -v unmappedpath="$1" "$AWK_SPLIT_UNMAPPED_SCRIPT" 
