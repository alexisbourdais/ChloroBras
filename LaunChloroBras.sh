#!  /usr/bin/bash

nextflow run ChloroBras.nf \
-profile slurm,singularity \
--singularity "-B /home:/home" \
--formatReadName .fq.gz \
--workflow fromReads \
--annotation all
#--qc \
#--quast \
#--trimming fastp \
#--seqtkSubsamp 0
