#!  /usr/bin/bash

. /local/env/envnextflow-23.10.0.sh

nextflow run ChloroBras.nf \
-profile standard,conda \
--formatReadName .fq.gz \
--workflow fromAsm \
--annotation organnot \
--singularity "-B /home:/home"
#--qc \
#--quast \
#--trimming fastp \
#--seqtkSubsamp 0 \
#--annotation all \