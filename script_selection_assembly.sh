#!/bin/bash

"""
This script allows to select the assembly of GetOrganelle whose ndhF gene is in the forward direction
The assemblies obtained for the same sample are specified in argument 1 and 2
Assembly for which the sequence specified in the pat variable is found is redirected to the file in argument 3
"""

#pattern in sequence of ndhF forward
pat="AAAGAAACTTGT"

#if pattern present, return 1 for the good assembly
res_1=$(grep -c ${pat} $1)
echo "$1 : $res_1"

res_2=$(grep -c ${pat} $2)
echo "$2 : $res_2"

#copy the good assembly in the output
if [[ $res_1 == 1 ]]; then
    while read a
    do
        echo $a
    done < $1 > $3
fi

if [[ $res_2 == 1 ]]; then
    while read a
    do
        echo $a
    done < $2 > $3
fi
