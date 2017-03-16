awk '/^#/ {print $0;next} {print $0 | "sort -k 2,2 -k 5,5 -k 3,3n -k 6,6n -k 1,1 --field-separator=\v"}'
