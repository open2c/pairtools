#!/usr/bin/env bash

INPUT='-'
OUTPUT='-'
OUTPUT_REST='-'
USAGE_LINE='Usage: pairsam_select_pair_type [-i input_file] [-o output_file] [-r output_rest_file] pair_type

where:
    -i  input file. By default, the input is read from stdin.
    -o  output file. By default, the output is printed into stdout.
    -r  output file for pairs of other types. By default, such pairs are dropped.
    pair_type  a quoted line specifying the requested pair types. Multiple pair types 
        should be separated by a | char. Examples: "LL", "CX|LL".
'

while getopts ":i:o:r:h" opt; do
  case $opt in
    i)
      INPUT=$OPTARG
      ;;
    o)
      OUTPUT=$OPTARG
      ;;
    r)
      OUTPUT_REST=$OPTARG
      ;;
    h)
      echo 'Read a pairsam file and print only the pairs of a certain type(s).'
      echo "${USAGE_LINE}" >&2
      exit 0
      ;;
    \?)
      echo "Invalid Option: -$OPTARG" 1>&2
      echo "${USAGE_LINE}" >&2
      exit 1
      ;;
  esac
done

shift $((OPTIND -1))

if [ $# -eq 1  ] ; then
    PAIRTYPE=$1
else
    echo "Please, provide the required pair type or types as the last argument!" 1>&2
    echo "${USAGE_LINE}"
    exit 1
fi

AWK_SCRIPT='
    BEGIN { FS="\v" } 
    {
        if ( $0 ~ /^#/ ) {print $0 '
if [[ "${OUTPUT}" != '-' ]]; then
    AWK_SCRIPT="${AWK_SCRIPT}"'> OUTPUT'
fi

if [[ "${OUTPUT_REST}" != '-' ]]; then
    AWK_SCRIPT="${AWK_SCRIPT}"'; print > OUTPUT_REST'
fi

AWK_SCRIPT="${AWK_SCRIPT}"' } else if ( $8 ~ /'"${PAIRTYPE}"'/ ) {print '

if [[ "${OUTPUT}" != '-' ]]; then
    AWK_SCRIPT="${AWK_SCRIPT}"'> OUTPUT'
fi

AWK_SCRIPT="${AWK_SCRIPT}"' } '

if [[ "${OUTPUT_REST}" != '-' ]]; then
    AWK_SCRIPT="${AWK_SCRIPT}"'else {print > OUTPUT_REST}'
fi

AWK_SCRIPT="${AWK_SCRIPT}"' } '

if [[ "${INPUT}" == '-' ]]; then
    awk -v OUTPUT=$OUTPUT -v OUTPUT_REST=$OUTPUT_REST "$AWK_SCRIPT" 
else
    awk -v OUTPUT=$OUTPUT -v OUTPUT_REST=$OUTPUT_REST "$AWK_SCRIPT" $INPUT
fi

