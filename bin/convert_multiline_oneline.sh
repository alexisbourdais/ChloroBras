#!/usr/bin/env bash

#Transforms a multi-line fasta into a single line.
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $1 > $2
