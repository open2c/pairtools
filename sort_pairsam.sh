[ $# -ge 1 -a -f "$1" ] && INPUT="$1" || INPUT="-"
cat $INPUT |\
    awk '/^#/ {print $0;next} {print $0 | "sort -k 1,1 -k 4,4 -k 2,2n -k 5,5n -k 8,8 --field-separator=\v"}'
