#!/bin/bash
if [ $# -ge 1 -a -f "$1" ] ; then
    awk '{if(/^#/)print;else exit}' $1
else 
    awk '{if(/^#/)print;else exit}'
fi 
