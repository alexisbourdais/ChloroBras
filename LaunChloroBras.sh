#!  /usr/bin/bash

. /local/env/envnextflow-23.10.0.sh

nextflow run ChloroBras.nf \
-profile slurm,conda \
--formatReadName .fq.gz \
--workflow assembling \
--assembler fastplast

#--singularity "-B /home:/home" \
#--qc \
#--quast \
#--trimming fastp \
#--seqtkSubsamp 0 \
#--annotation all \
#--phylogeny all
