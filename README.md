# ChloroBras!

## Overview

![Alt text](https://user-images.githubusercontent.com/108393526/176895796-71946738-161f-4c90-a7ec-bf478ae8bbcf.png)

ChloroBras is a project allowing the automatic assembly and analysis of chloroplast genome, developed for *Brassica* but transposable to any family of flowering plants. It consists of two nextflow pipelines.
The first one (**test_assembler.nf**) allows the assembly of chloroplastic genomes from paired Illumina reads, via three different assemblers (**GetOrganelle**, **Fast-Plast**, **ORGanelle ASseMbler**). It then aligns with **Nucmer** these genomes with a reference genome (**brassica_oleracea.fasta**) and allows the visualization of the quality of these assemblies via a dot-plot created by **Mummer**. A sub-sampling step via **Seqtk** was added for the **Fast-Plast** and **ORGanelle ASseMbler** assemblers, as the samples from which the pipeline was developed were originally intended for the study of nuclear polymorphisms, the assembly could take several days because of the large number of reads present. (**GetOrganelle** is able to perform its own subsampling.)

The second pipeline (**analysis_chloro.nf**) allows the assembly by **Getorganelle** (which showed the best results on our data), the selection by a bash script (**script_selection_assembly.sh**) of the assembly with the Small Single Copy in the right direction (**GetOrganelle** provides two assemblies per sample with the only difference being the direction of the SSC), the alignment of the different assemblies by **Mafft** and finally the creation of a phylogenetic tree by **RAxML**.

A python script (**rename_fasta_header.py**) allows to rename the headers after the assembly according to the name of the sample and the assembler.
A bash script (**rename_fasta_header_phylo.sh**) that reads a csv file (**tab_correspondence_phylo.csv**) allows to rename the headers of the sequences in order to display the desired information on the phylogenetic tree.

## Installation

Nextflow and most tools can be installed from conda

### Nextflow

`conda install -c bioconda nextflow`

### Singularity

`conda install -c conda-forge singularity`

### Seqtk

`conda install -c bioconda seqtk`

### GetOrganelle

`conda install -c bioconda getorganelle`

### Mummer/Nucmer

`conda install -c bioconda mummer`

### Mafft

`conda install -c bioconda mafft`

### RAxML

`conda install -c bioconda mummer`

### Fast-Plast and ORGanelle ASseMbler

These tools have been integrated by a docker image loaded by singularity.

## Instruction

After loading the nextflow and singularity environment, you just have to use the following command to start the pipeline:
`nextflow run *pipeline*.nf`


### Paramaters


Defines the pipeline inputs parameters

Each of the following parameters can be specified as command line options or in the config file (**nextflow.config**)

Beware of quotes and do not forget to modify the path when necessary


##### Test_assembler

    Reads directory
    --readsFiles                    Path to input data, default: "./samples/*_R{1,2}.fastq.gz"

    Results directory
    --resultsDir                    Path to results directory, default: "./results/"

    Sqtk
    --seqtk_nb_read                 Number of reads to keep, default: 2000000

    GetOrganelle
    --getorganelle_index            Index of GetOrganelle, default: "embplant_pt"
    --getorganelle_kmer             Size of kmers, default: "21,45,65,85,105"

    FastPlast
    --fastplast_index               Index of Fast-Plast, default: "Brassicales"

    OrgAsm
    --orgasm_probes                 Index of ORGanelle ASeMbler, default: "protChloroArabidopsis"

    Rename_headers
    --rename_dir_script             Absolute path to script rename

    Nucmer
    --nucmer_ref                    Absolute path to Fasta reference for alignment, default: brassica oleracea

    Mummer
    --mummer_axe                    Size of X-axis (fonction of genome's size), default (plastome): "'[0:154000]'"
    --mummer_format_output          Format of the plot, default: "png"

##### Analysis_chloro

    Reads directory
    --readsFiles                    Path to input data, default: "./samples/*_R{1,2}.fastq.gz"

    Results directory
    --resultsDir                    Path to results directory, default: "./results/"

    GetOrganelle
    --getorganelle_index            Index of GetOrganelle, default: "embplant_pt"
    --getorganelle_kmer             Size of kmers, default: "21,45,65,85,105"

    Rename_headers
    --rename_dir_script             Absolute path to script rename

    Nucmer
    --nucmer_ref                    Absolute path to Fasta reference for alignment

    Mummer
    --mummer_axe                    Size of X-axis (fonction of genome's size), default (plastome): "'[0:154000]'"
    --mummer_format_output          Format of the plot, default: "png"

    Select_assembly
    --select_assembly_script        Absolute path to script select_assembly

    Mafft
    --mafft_method                  Alignment methods, default: "auto"

    Rename_headers_phylo
    --rename_headers_phylo_script   Absolute path to script rename

    Raxml
    --raxml_model                   Model uses by RAxML, default: "GTRGAMMAI"
    
