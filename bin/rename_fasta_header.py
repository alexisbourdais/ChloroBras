#!/usr/bin/env python3

'''
Change the header of the sequence contained in the file specified by the -i argument
The name of the desired header is to be specified with the -n argument 
The output file is specified by the -o argument
'''


import argparse
# Construct the argument parser
ap = argparse.ArgumentParser()
# Add the arguments to the parser
ap.add_argument("-i", "--infile", required=True, help="fasta file")
ap.add_argument("-n", "--name", required=True, help="name")
ap.add_argument("-o", "--output", required=True, help="fasta file output")

args = vars(ap.parse_args())
seq=0
with open(args['infile'], "r") as fasta_in:
    with open(args['output'], "w+") as fasta_out:
        for line in fasta_in:
            if line[0] == ">":
                #seq=seq+1
                fasta_out.write('>'+args["name"]+'\n')
            else:
                fasta_out.write(line)


#+"_"+str(seq)+
