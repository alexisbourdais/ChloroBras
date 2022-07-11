# ChloroBras!

## Overview

![Alt text](https://user-images.githubusercontent.com/108393526/176895796-71946738-161f-4c90-a7ec-bf478ae8bbcf.png)

ChloroBras is a project allowing the automatic assembly and analysis of chloroplast genome, developed for *Brassica* but transposable to any family of flowering plants. It consists of two nextflow pipelines and several bash or python scripts.

The first one (**test_assembler.nf**) allows the assembly of chloroplastic genomes from paired Illumina reads, via three different assemblers (**GetOrganelle**, **Fast-Plast**, **ORGanelle ASseMbler**). It then aligns with **Nucmer** these genomes with a reference genome (**brassica_oleracea.fasta**) and allows the visualization of the quality of these assemblies via a dot-plot created by **Mummer**. A sub-sampling step via **Seqtk** was added for the **Fast-Plast** and **ORGanelle ASseMbler** assemblers, as the samples from which the pipeline was developed were originally intended for the study of nuclear polymorphisms, the assembly could take several days because of the large number of reads present. (**GetOrganelle** is able to perform its own subsampling.)

The second pipeline (**analysis_chloro.nf**) allows the assembly by **Getorganelle** (which showed the best results on our data), the selection by a bash script (**script_selection_assembly.sh**) of the assembly with the Small Single Copy in the right direction (**GetOrganelle** provides two assemblies per sample with the only difference being the direction of the SSC), the alignment of the different assemblies by **Mafft** and finally the creation of a phylogenetic tree by **RAxML**.

The bash script (**script_selection_assembly.sh**) selects the correctly structured GetOrganelle assembly thanks to a short highly conserved sequence of the ndhF gene located on the SSC.

The python script (**rename_fasta_header.py**) renames the headers after the assembly according to the name of the sample and the assembler.

The bash script (**rename_fasta_header_phylo.sh**) that reads a csv file (**tab_correspondence_phylo.csv**) allows to rename the headers of the sequences in order to display the desired information on the phylogenetic tree.


## Installation

Nextflow and singularity can be installed from conda. All tool environments using conda in the pipeline will be automatically created and activated.

### Nextflow

`conda install -c bioconda nextflow`

### Singularity

`conda install -c conda-forge singularity`

### Fast-Plast and ORGanelle ASseMbler

These tools have been integrated by a docker image loaded by singularity. (https://hub.docker.com/u/chloroextractorteam)

## Instruction

Before starting a pipeline, a **directory with all the samples must be prepared**. These should have a structured name like: xxx_R1.fastq.gz and xxx_R2.fastq.gz. It is possible to use symbolic links, which can be created with the following command:
`ln -s path/to/xxx_R1.fastq.gz xxx_R1.fastq.gz`

The **csv file must also be prepared in advance**. The csv file must also be prepared in advance for the analysis_chloro pipeline. It is composed of the name (xxx) of the sample in the first column, the name of the assemblies obtained by GetOrganelle in the 2nd (xxx_get_1_1) and 3rd column (xxx_get_2_1), and the name to display in the 4th column.


After loading the nextflow and singularity environment, you just have to use the following command to start the pipeline:
`nextflow run *pipeline*.nf --option` (Options are optional, see next topic) 


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
    
## Documentation

Fast-Plast: https://github.com/mrmckain/Fast-Plast

GetOrganelle: https://github.com/Kinggerm/GetOrganelle

MAFFT: https://mafft.cbrc.jp/alignment/software/

MUMmer/NUCMER: https://mummer4.github.io/index.html

Nextflow: https://www.nextflow.io/docs/latest/index.html

ORGanelle ASseMbler: https://git.metabarcoding.org/org-asm/org-asm

RAxML: https://cme.h-its.org/exelixis/web/software/raxml/

Seqtk: https://github.com/lh3/seqt

## References
    
P. Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017) doi:10.1038/nbt.3820

Freudenthal, Jan A., Simon Pfaff, Niklas Terhoeven, Arthur Korte, Markus J. Ankenbrand, et Frank Förster. « A Systematic Comparison of Chloroplast Genome Assembly Tools ». Genome Biology 21, n o 1 (décembre 2020): 254. https://doi.org/10.1186/s13059-020-02153-6.

Jin, Jian-Jun, Wen-Bin Yu, Jun-Bo Yang, Yu Song, Claude W. dePamphilis, Ting-Shuang Yi, et De-Zhu Li.« GetOrganelle: A Fast and Versatile Toolkit for Accurate de Novo Assembly of Organelle Genomes ». GenomeBiology 21, n o 1 (décembre 2020): 241. https://doi.org/10.1186/s13059-020-02154-5.

Katoh, K. « MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform ». Nucleic Acids Research 30, n o 14 (15 juillet 2002): 3059-66. https://doi.org/10.1093/nar/gkf436.

Kurtzer, Gregory M., Vanessa Sochat, et Michael W. Bauer. « Singularity: Scientific Containers for Mobility of Compute ». Édité par Attila Gursoy. PLOS ONE 12, n o 5 (11 mai 2017): e0177459. https://doi.org/10.1371/journal.pone.0177459.

Marçais, Guillaume, Arthur L. Delcher, Adam M. Phillippy, Rachel Coston, Steven L. Salzberg, et Aleksey Zimin. « MUMmer4: A Fast and Versatile Genome Alignment System ». Édité par Aaron E. Darling. PLOS Computational Biology 14, n o 1 (26 janvier 2018): e1005944. https://doi.org/10.1371/journal.pcbi.1005944.

Stamatakis, Alexandros. « RAxML Version 8: A Tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies ». Bioinformatics 30, n o 9 (1 mai 2014): 1312-13. https://doi.org/10.1093/bioinformatics/btu033.

Fast-Plast: McKain et Wilson, 2017

ORGanelle ASseMbler : Coissac et al. 2019
