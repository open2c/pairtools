#!/bin/bash
if [ $# -ge 1 -a -f "$1" ] ; then
    sed -n -e '/^[^#]/,$p' $1
else 
    sed -n -e '/^[^#]/,$p'
fi 

