#!  /usr/bin/bash

. /local/env/envnextflow-23.10.0.sh

nextflow run ChloroBras.nf \
--formatReadName .fq.gz \
--workflow assembling \
-profile slurm,singularity \
--assembler getorganelle