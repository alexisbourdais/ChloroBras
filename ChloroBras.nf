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

    Command : nextflow run ChloroBras.nf -profile [standard/slurm,singularity/conda] --workflow [fromReads/fromAsm] --singularity "-B root/to/mount/"

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
    --refFasta              Path to Fasta reference for alignment and quast, default: "./Data/Brassica-oleracea-isolate-HDEM-chloroplast.fasta"
    --refGff                Path to Gff reference for quast, default: "./Data/Brassica-oleracea-isolate-HDEM-chloroplast.gff3"
    --mummerAxe             Size of X-axis (fonction of genome's size), default (plastome): "'[0:154000]'"
    --mummerFormatOut       Format of the plot, default: "png"

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

///////////////////////////////////////////////////////////
//////////////////////     Process      ///////////////////
///////////////////////////////////////////////////////////

results = file(params.resultsDir)

//////////////////////  QC  //////////////////////

process fastqc {

    label 'process_low'

    input:
    path(read)

    output :
    path("${read.simpleName}_fastqc.zip")

    script:
    """
    fastqc -t ${task.cpus} ${read} -o "./"
    """
}

process multiqc {

    label 'process_low'

    publishDir "${results}/", mode: 'copy'

    input:
    path(allFastq)

    output :
    path("multiqc_report.html")

    script:
    """
    multiqc ${allFastq}
    """
}

//////////////////////  TRIMMING  //////////////////////

process trimgalore {

    label 'process_low'

    tag "${sampleId}"

    publishDir "${results}/Trimming_trimgalore/", mode: 'copy'

    input:
    tuple val(sampleId), path(reads)

    output :
    tuple val("${sampleId}"), path("${sampleId}_R{1,2}_trimgalore.fq.gz"), emit: paired_reads
    tuple path("${sampleId}_R1_trimgalore_fastqc.zip"), path("${sampleId}_R2_trimgalore_fastqc.zip"), emit: qc
    tuple path("${sampleId}_R1_trimgalore_report.txt"), path("${sampleId}_R2_trimgalore_report.txt")

    script:
    """
    trim_galore \
    --paired ${reads[0]} ${reads[1]} \
    --basename ${sampleId} \
    --gzip \
    -o "Trimming_trimgalore" \
    --cores "${task.cpus}" \
    --fastqc
    #--illumina --stranded_illumina

    mv "Trimming_trimgalore/${sampleId}_val_1.fq.gz" "${sampleId}_R1_trimgalore.fq.gz"
    mv "Trimming_trimgalore/${sampleId}_val_2.fq.gz" "${sampleId}_R2_trimgalore.fq.gz"
    mv "Trimming_trimgalore/${sampleId}_val_1_fastqc.zip" "${sampleId}_R1_trimgalore_fastqc.zip"
    mv "Trimming_trimgalore/${sampleId}_val_2_fastqc.zip" "${sampleId}_R2_trimgalore_fastqc.zip"
    mv "Trimming_trimgalore/${reads[0]}_trimming_report.txt" "${sampleId}_R1_trimgalore_report.txt"
    mv "Trimming_trimgalore/${reads[1]}_trimming_report.txt" "${sampleId}_R2_trimgalore_report.txt"
    """
}

process fastp {

    label 'process_low'

    tag "${sampleId}"

    publishDir "${results}/Trimming_fastp/", mode: 'copy'

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val("${sampleId}"), path("${sampleId}_R{1,2}_trimfastp.fq.gz"), emit: paired_reads
    path("report_fastp.html")

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
    -o ${sampleId}_R1_trimfastp.fq.gz -O ${sampleId}_R2_trimfastp.fq.gz \
    --thread ${task.cpus} \
    --html report_fastp.html
    """
}

//////////////////////  ASSEMBLER  //////////////////////

process getorganelle_index {

    input:
    val index

    output:
    val("")

    script:
    """
    get_organelle_config.py -a ${index}
    """
}

process getorganelle {

    label 'process_high'

    tag "${sampleId}"

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    publishDir "${results}/Assembly/", mode: 'copy', pattern: "*getGood.fasta"

    input:
    val index
    tuple val(sampleId), path(reads)

    output :
    path("${sampleId}_getGood.fasta"), emit: mafft
    tuple val(sampleId), val("get"), path("${sampleId}_get.fasta"), emit: mummer

    script:
    """
    get_organelle_from_reads.py \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -o ${sampleId} \
    -k ${params.getKmer} \
    -F ${params.getIndex} \
    -t ${task.cpus}

    asm1="\$(ls ${sampleId}/embplant_pt.K*.complete.graph1.1.path_sequence.fasta)"
    asm2="\$(ls ${sampleId}/embplant_pt.K*.complete.graph1.2.path_sequence.fasta)"

    rename_fasta_header.py -i "\${asm1}" -n "${sampleId}_get_1" -o "${sampleId}_get_1.fasta"
    rename_fasta_header.py -i "\${asm2}" -n "${sampleId}_get_2" -o "${sampleId}_get_2.fasta"
    cat "${sampleId}_get_1.fasta" "${sampleId}_get_2.fasta" > "${sampleId}_get.fasta"
    selection_assembly.sh "${sampleId}_get_1.fasta" "${sampleId}_get_2.fasta" "${sampleId}_getGood.fasta"
    """
}

process seqtk {

    tag "${sampleId}"
    
    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), path("${sampleId}_sub{1,2}.fq.gz")

    script:
    """
    seqtk sample -s666 ${reads[0]} ${params.seqtkSubsamp} |gzip -c - > ${sampleId}_sub1.fq.gz
    seqtk sample -s666 ${reads[1]} ${params.seqtkSubsamp} |gzip -c - > ${sampleId}_sub2.fq.gz
    """
}

process fastplast {

    label 'process_high'

    publishDir "${results}/Assembly/", mode: 'copy'

    tag "${sampleId}"

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), val("fpt"), path("${sampleId}_fpt.fasta"), emit: mummer
    path("${sampleId}_fpt.fasta"), emit: mafft

    script:
    """
    fast-plast.pl \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    --name ${sampleId} \
    --bowtie_index ${params.fastIndex} \
    --threads ${task.cpus}

    mv "${sampleId}/Final_Assembly/${sampleId}_FULLCP.fsa" "${sampleId}_fpt.fasta"
    """
}

process orgasm {

    label 'process_high'

    publishDir "${results}/Assembly/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    tag "${sampleId}"

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), val("org"), path("${sampleId}_org.fasta"), emit: mummer
    path("${sampleId}_org.fasta"), emit: mafft

    script:
    """
    oa index --estimate-length=0.9 ${sampleId} ${reads[0]} ${reads[1]}
    oa buildgraph --probes ${params.orgasmProbes} ${sampleId} ${sampleId}
    oa unfold ${sampleId} ${sampleId} > ${sampleId}_orgasm.fasta
    convert_multiline_oneline.sh "${sampleId}_orgasm.fasta" "${sampleId}_oneLine_org.fasta"
    rename_fasta_header.py -i "${sampleId}_oneLine_org.fasta" -n "${sampleId}" -o "${sampleId}_org.fasta"
    """
}

//////////////////////  QUALITY  //////////////////////

process quast {

    label 'process_low'

    publishDir "${results}/", mode: 'copy'

    input:
    val(asm)
    path(allAsm)

    output:
    tuple path("Quast_${asm}/transposed_report.tsv"), path("Quast_${asm}/report.html"), path("Quast_${asm}/report.html"), path("Quast_${asm}/circos/circos.png"), path("Quast_${asm}/genome_stats/"), path("Quast_${asm}/predicted_genes/*.gff")

    script:
    """
    quast -o Quast_${asm} -t ${task.cpus} -r ${params.refFasta} --circos --glimmer --features ${params.refGff} ${allAsm}
    """
}

process mummer {

    label 'process_low'

    publishDir "${results}/Mummer", mode: 'copy'

    input:
    tuple val(sampleId), val(assembler), path(assembly)

    output:
    path("${sampleId}_${assembler}.png")

    script:
    """
    nucmer -p ${sampleId} ${params.refFasta} ${assembly}

    mummerplot \
    -x ${params.mummerAxe} \
    -p ${sampleId}_${assembler} \
    -t ${params.mummerFormatOut} \
    "${sampleId}.delta"
    """
}

//////////////////////  ANNOTATION  //////////////////////

process mfannot {

    label 'process_medium'

    cache false

    publishDir "${results}/Annotation/", mode: 'move'

    input:
    tuple val(sampleId), val(assembler), path(assembly)

    output:
    path("${sampleId}_${assembler}_mfannot.txt")

    script:
    """
    mfannot -g 11 -o ${sampleId}_${assembler}_mfannot.txt ${assembly}
    """
}

process organnot {

    label 'process_medium'

    publishDir "${results}/Annotation/", mode: 'copy'

    input:
    tuple val(sampleId), val(assembler), path(assembly)

    output:
    path("${sampleId}_${assembler}_organnot.txt")

    script:
    if (params.workflow == "fromAsm") {
        """
        organnot -c ${assembly}
        mv "${sampleId}.annot.circular" "${sampleId}_${assembler}_organnot.txt"
        """
    }
    else {
        """
        organnot -c ${assembly}
        mv ${sampleId}_${assembler}.annot.circular ${sampleId}_${assembler}_organnot.txt
        """
    }

}

//////////////////////  ALIGNER  //////////////////////

process mafft {

    label 'process_high'

    publishDir "${results}/Mafft/", mode: 'copy'

    input:
    path(multiFasta)

    output:
    path("align.fasta")

    script:
    """
    mafft --${params.mafftMethod} ${multiFasta} > align.fasta
    sed -i -e 's/_get_1//g' -e 's/_get_2//g' align.fasta
    """
}

//////////////////////  PHYLOGENY  //////////////////////

process raxml {

    label 'process_high'

    errorStrategy 'ignore'

    publishDir "${results}/RAxML/", mode: 'copy'

    input:
    path(multi_fasta_align)

    output:
    path("*.model_${params.raxmlModel}")

    script:
    """
    raxmlHPC \
    -s ${multi_fasta_align} \
    -n model_${params.raxmlModel} \
    -d \
    -m ${params.raxmlModel} \
    --auto-prot=bic \
    --bootstop-perms=1000 \
    -x 12345 \
    -p 12345 \
    -# autoMRE \
    -f a
    """
}

process iqtree {

    label 'process_high'

    errorStrategy 'ignore'

    publishDir "${results}/Iqtree/", mode: 'copy'

    input:
    path(multi_fasta_align)

    output:
    path("iqtree_${params.iqtreeModel}.*")

    script:
    if (params.iqtreeOption=="") {
    """
    iqtree \
    -s ${multi_fasta_align} \
    -pre "iqtree_${params.iqtreeModel}" \
    -m ${params.iqtreeModel} \
    -T ${task.cpus}
    """
    }

    else {
    """
    iqtree \
    -s ${multi_fasta_align} \
    -pre "iqtree_${params.iqtreeModel}" \
    -m ${params.iqtreeModel} \
    -T ${task.cpus} \
    ${params.iqtreeOption}
    """
    }
}

process raxmlng {

    label 'process_high'

    errorStrategy 'ignore'

    publishDir "${results}/RAxML-NG/", mode: 'copy'

    input:
    path(multi_fasta_align)

    output:
    path("raxmlng_${params.raxmlngModel}.raxml.*")

    script:

    if (params.raxmlngOption=="") {
    """
    raxml-ng \
    --all \
    --bs-trees ${params.raxmlngBootstrap} \
    --msa ${multi_fasta_align} \
    --model ${params.raxmlngModel} \
    --prefix "raxmlng_${params.raxmlngModel}" \
    --threads ${task.cpus}
    """
    }

    else {
    """
    raxml-ng \
    --all \
    --bs-trees ${params.raxmlngBootstrap} \
    --msa ${multi_fasta_align} \
    --model ${params.raxmlngModel} \
    --prefix "raxmlng_${params.raxmlngModel}" \
    --threads ${task.cpus} \
    ${params.raxmlngOption==""}
    """
    }
}

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

workflow.onComplete{println("Workflow execution completed sucessfully ! or not ...")}
