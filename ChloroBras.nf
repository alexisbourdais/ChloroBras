#!/local/env/env nextflow

/*
===============================================================
 Plastome Analysis Pipeline. Started April 2022.
 #### Homepage / Documentation
 https://github.com/alexisbourdais/
 #### Authors
 Alexis Bourdais
---------------------------------------------------------------
*/

///////////////////////////////////////////////////////////
//////////////////////       HELP       ///////////////////
///////////////////////////////////////////////////////////

def helpMessage() {
    log.info """

    Command : nextflow run ChloroBras.nf -profile [standard/slurm,singularity/conda] --workflow [fromReads/fromAsm] --singularity "-B /home:/home"

    REQUIRED parameter

    -profile [standard]/slurm,      Select profile standard (local) or slurm. Default: standard          
             singularity/conda      Select profile singularity or conda. (FastPlast, Orgasm, mfannot and organnot are only available with singularity, even in conda profile)
                                                                         (Mummer is only available with conda, even in singularity profile)

    --workflow [fromReads/fromAsm]     fromReads : chloroplast genome assembly, annotation, quality assessment with quast and dotplot, phylogeny analysis from paired reads
                                       fromAsm : mafft alignement, annotation, phylogeny analysis from assemblies

    Singularity
    --singularity           Mounted directory, default: "-B /scratch:/scratch -B /home:/home -B /local:/local -B /db:/db -B /groups:/groups"

    OPTIONAL parameter              

    Reads directory
    --readDir               Default: "./Data"
    --baseReadName          Default: "_R{1,2}"     ex: name_R1.fastq.gz & name_R2.fastq.gz
    --formatReadName        Default: ".fastq.gz"

    Results directory
    --resultsDir            Path to results directory, default: "./Results/"

    Assembly directory
    --assemblyDir           Path to assembly directory, default: "./Results/Assembly/"
    --formatAsm             Default: ".fasta"

    Assembler
    --assembler             Choose assembler to use (['getorganelle'], 'fastplast', 'orgasm' or 'all') 
                            Phylogeny analysis is not available with 'orgasm' and 'all'.
    Quality control
    --qc                    To activate qc

    Trimming
    --trimming              Add trimming step with 'fastp' or 'trimgalore'. Default: none

    Annotation
    --annotation            Add annotation step ('all', 'mfannot', 'organnot'). Default: none            

    GetOrganelle
    --getIndex              Index of GetOrganelle, default: "embplant_pt"
    --getKmer               Size of kmers, default: "21,45,65,85,105"

    Sqtk
    --seqtkSubsamp          Subsampling for Orgasm and FastPlast, default: 2000000. 
                            Set to 0 to deactivate (assembly can be time-consuming)
    FastPlast
    --fastIndex             Index of Fast-Plast, default: "Brassicales"

    OrgAsm
    --orgasmProbes          Index of ORGanelle ASeMbler, default: "protChloroArabidopsis"

    Mummer - Quast
    --quast                 Activate quast : produce stats and circos between ref and assemblies.
    --refGff                Path to Gff reference for quast, default: "./Data/Brassica-oleracea-isolate-HDEM-chloroplast.gff3"
    --refFasta              Path to Fasta reference for alignment and quast, default: "./Data/Brassica-oleracea-isolate-HDEM-chloroplast.fasta"

    Mafft
    --mafftMethod           Alignment methods, default: "auto"

    Phylogeny
    --phylogeny             Add phylogenetic step ('raxml', 'iqtree', 'raxmlng' or 'all'). Default : none

    Raxml
    --raxmlModel            Model uses by RAxML, default: "GTRGAMMAI"

    IQtree
    --iqtreeModel           Model uses by IQtree, default: "GTR+I+G"
    --iqtreeOption          Use to add option to iqtree: "--option argument"

    Raxml-ng
    --raxmlngModel          Model uses by RAxML-NG, default: "GTR+G+I"
    --raxmlngBootstrap      Bootstrap number, default: 200
    --raxmlngOption         Use to add option to Raxml-ng: "--option argument"

    Each of the previous parameters can be specified as command line options, in launch file or in the config file

    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
        helpMessage()
        exit 0
}

//////////////////////////////////////////////
////////////    Modules     //////////////////
//////////////////////////////////////////////

include { fastp }                               from './modules/fastp'
include { trimgalore }                          from './modules/trimgalore'
include { fastqc }                              from './modules/fastqc'
include { multiqc }                             from './modules/multiqc'
include { seqtk }                               from './modules/seqtk'
include { getorganelle_index; getorganelle }    from './modules/getorganelle'
include { fastplast }                           from './modules/fastplast'
include { orgasm }                              from './modules/orgasm'
include { quast }                               from './modules/quast'
include { mummer }                              from './modules/mummer'
include { mfannot }                             from './modules/mfannot'
include { organnot }                            from './modules/organnot'
include { mafft }                               from './modules/mafft'
include { raxml }                               from './modules/raxml'
include { raxmlng }                             from './modules/raxmlng'
include { iqtree }                              from './modules/iqtree'

///////////////////////////////////////////////////////////
//////////////////////     Sub-Workflow     ///////////////
///////////////////////////////////////////////////////////

workflow quality_wf {
    take:
    reads

    main:
    fastqc(reads)
    multiqc(fastqc.out.collect())
}

workflow getorganelle_wf {
    take:
    paired_reads

    main:
    getorganelle_index(params.getIndex)
    getorganelle(getorganelle_index.out, paired_reads)
    mummer(getorganelle.out.mummer)

    if (params.quast) {
        quast("get", getorganelle.out.mafft.collect())
    }

    if (params.annotation=="mfannot" || params.annotation=="all") {
        mfannot(getorganelle.out.mummer)
    }
    
    if (params.annotation=="organnot" || params.annotation=="all") {
        organnot(getorganelle.out.mummer)
    }

    emit:
    assembly = getorganelle.out.mafft
}

workflow fastplast_wf {
    take:
    sub_paired_reads

    main:
    fastplast(sub_paired_reads)
    mummer(fastplast.out.mummer)

    if (params.quast) {
        quast("fpt", fastplast.out.mafft.collect())
    }

    if (params.annotation=="mfannot" || params.annotation=="all") {
        mfannot(fastplast.out.mummer)
    }
    
    if (params.annotation=="organnot" || params.annotation=="all") {
        organnot(fastplast.out.mummer)
    }

    emit:
    assembly = fastplast.out.mafft
}

workflow orgasm_wf {
    take:
    sub_paired_reads

    main:
    orgasm(sub_paired_reads)
    mummer(orgasm.out.mummer)

    if (params.quast) {
        quast("org", orgasm.out.mafft.collect())
    }

    if (params.annotation=="mfannot" || params.annotation=="all") {
        mfannot(orgasm.out.mummer)
    }

    if (params.annotation=="organnot" || params.annotation=="all") {
        organnot(orgasm.out.mummer)
    }

    emit:
    assembly = orgasm.out.mafft
}

workflow phylogeny_wf {
    take:
    assembly

    main:
    mafft(assembly.collectFile(name: 'multi_fasta', newLine: true))

    if (params.phylogeny=="raxml" || params.phylogeny=="all") {
        raxml(mafft.out)
    }

    if (params.phylogeny=="iqtree" || params.phylogeny=="all") {
        iqtree(mafft.out)
    }

    if (params.phylogeny=="raxmlng" || params.phylogeny=="all") {
        raxmlng(mafft.out)
    }
}

///////////////////////////////////////////////////////
//////////////////     Workflow     ///////////////////
///////////////////////////////////////////////////////


workflow {

    if (params.workflow=="fromReads") {

        paired_reads = Channel.fromFilePairs("${params.readFiles}", checkIfExists:true)
        reads = Channel.fromPath("${params.readFiles}")

        if (params.qc) {
            quality_wf(reads)
        }

        // GetOrganelle
        if (params.assembler=='getorganelle') {

            // Trimming Fastp
            if (params.trimming=="fastp") {
                fastp(paired_reads)
                getorganelle_wf(fastp.out.paired_reads)
            }

            // Trimming Trimgalore
            else if (params.trimming=="trimgalore") {
                trimgalore(paired_reads)
                getorganelle_wf(trimgalore.out.paired_reads)
            }

            // Without trimming
            else {
                getorganelle_wf(paired_reads) 
            }

            if (params.phylogeny) {
                phylogeny_wf(getorganelle_wf.out.assembly)
            }
        }

        // FastPlast
        else if (params.assembler=='fastplast') {

            // Trimming Fastp
            if (params.trimming=="fastp") {
                fastp(paired_reads)

                // Without subsampling
                if (params.seqtkSubsamp==0) {
                    fastplast_wf(fastp.out.paired_reads)
                }

                // With subsampling
                else {
                    seqtk(fastp.out.paired_reads)
                    fastplast_wf(seqtk.out)
                }
            }

            // Trimming Trimgalore
            else if (params.trimming=="trimgalore") {
                trimgalore(paired_reads)

                // Without subsampling
                if (params.seqtkSubsamp==0) {
                    fastplast_wf(trimgalore.out.paired_reads)
                }

                // With subsampling
                else {
                    seqtk(trimgalore.out.paired_reads)
                    fastplast_wf(seqtk.out)
                }
            }

            // Without trimming
            else {

                // Without subsampling
                if (params.seqtkSubsamp==0) {
                    fastplast_wf(paired_reads)
                }

                // With subsampling
                else {
                    seqtk(paired_reads)
                    fastplast_wf(seqtk.out)
                }
            }

            if (params.phylogeny) {
                phylogeny_wf(fastplast_wf.out.assembly)
            }
        }

        // Orgasm
        else if (params.assembler=='orgasm') {

            // Trimming Fastp
            if (params.trimming=="fastp") {
                fastp(paired_reads)

                // Without subsampling
                if (params.seqtkSubsamp==0) {
                    orgasm_wf(fastp.out.paired_reads)
                }

                // With subsampling
                else {
                    seqtk(fastp.out.paired_reads)
                    orgasm_wf(seqtk.out)
                }
            }

            // Trimming Trimgalore
            else if (params.trimming=="trimgalore") {
                trimgalore(paired_reads)

                // Without subsampling
                if (params.seqtkSubsamp==0) {
                    orgasm_wf(trimgalore.out.paired_reads)
                }

                // With subsampling
                else {
                    seqtk(trimgalore.out.paired_reads)
                    orgasm_wf(seqtk.out)
                }
            }

            // Without trimming
            else {
                // Without subsampling
                if (params.seqtkSubsamp==0) {
                    orgasm_wf(paired_reads)
                }

                // With subsampling
                else {
                    seqtk(paired_reads)
                    orgasm_wf(seqtk.out)
                }
            }
        }

        // All
        else if (params.assembler=='all') {

            // Trimming Fastp
            if (params.trimming=="fastp") {
                fastp(paired_reads)
                getorganelle_wf(fastp.out.paired_reads)

                // Without subsampling
                if (params.seqtkSubsamp==0) {
                    fastplast_wf(fastp.out.paired_reads)
                    orgasm_wf(fastp.out.paired_reads)
                }

                // With subsampling
                else {
                    seqtk(fastp.out.paired_reads)
                    fastplast_wf(seqtk.out)
                    orgasm_wf(seqtk.out)
                }
            }

            // Trimming Trimgalore
            else if (params.trimming=="trimgalore") {
                trimgalore(paired_reads)
                getorganelle_wf(trimgalore.out.paired_reads)

                // Without subsampling
                if (params.seqtkSubsamp==0) {
                    fastplast_wf(trimgalore.out.paired_reads)
                    orgasm_wf(trimgalore.out.paired_reads)
                }

                // With subsampling
                else {
                    seqtk(trimgalore.out.paired_reads)
                    fastplast_wf(seqtk.out)
                    orgasm_wf(seqtk.out)
                }
            }

            // Without trimming
            else {
                getorganelle_wf(paired_reads)

                // Without subsampling
                if (params.seqtkSubsamp==0) {
                    fastplast_wf(paired_reads)
                    orgasm_wf(paired_reads)
                }

                // With subsampling
                else {
                    seqtk(paired_reads)
                    fastplast_wf(seqtk.out)
                    orgasm_wf(seqtk.out)
                }
            }
        } 
    }

    else if (params.workflow=="fromAsm") { 
        assemblies = Channel.fromPath("${params.asmFiles}")

        if (params.phylogeny) {
            phylogeny_wf(assemblies)
        }

        if (params.annotation) {
            name = assemblies.map { it.baseName }
            assembler = Channel.of("preAsm")
            input = name.merge(assembler).merge(assemblies)

            if (params.annotation=="mfannot" || params.annotation=="all") {
                mfannot(input)
            }
            if (params.annotation=="organnot" || params.annotation=="all") {
                organnot(input)
            }
        }
    }

    else { 
        println """
    
        Error : Any workflow selected !
    
        """
    helpMessage()
    exit 0 
    }
}

workflow.onComplete{
    if (workflow.success) {
        println("Workflow execution completed sucessfully !")
    } 
    else {
        println("Workflow execution completed with errors !")
    }
}
