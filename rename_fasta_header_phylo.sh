#!/bin/bash

"""
This script allows to modify the headers of a multi-fasta file specified as first argument
The current and desired headers must be saved in a vcf file
The corresponding headers must be on the same line, with the desired header in the 4th field
All sequences with their new headers will be saved in an output file specified as second argument
"""

ficcsv=./tab_correspondence_phylo.csv


#####################
#    Main program   #
#####################

ficin=$1
pat='^>.*'

while read a
do
	if [[ $a =~ $pat ]]; then
		rech=$(echo $a | sed "s/>//")
		sub=$(grep -w $rech $ficcsv | cut -d\; -f4)
		if [ -z "$sub" ]; then
			echo $a
		else
			echo ">$sub"
		fi
	else
		echo $a
	fi
done < $ficin > $2
