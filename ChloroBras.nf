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

    Command : nextflow run ChloroBras.nf --workflow [test/analysis/fromAsm]

    REQUIRED parameter

    Workflow
    --workflow [assembling/analyzing/fromAsm]       assembling : assembles genomes with the three assemblers, allows quality assessment via dotplot
                                                    analyzing : assemble genomes with [getorganelle] or Fastplast and create phylogenetic tree
                                                    fromAsm : mafft alignement and Raxml tree from assemblies in ./Results/Assembly/
    
    OPTIONAL parameter

    Executor
    --executor              Choose executor (local or slurm). Default: local

    Assembler
    --assembler             Choose assembler to use (getorganelle, fastplast or orgasm), default: all for assembling workflow
                                                                                         default: getorganelle for analysing workflow
    
    Reads directory
    --readDir                Default: "./Data"
    --baseReadName           Default: "_R{1,2}"
    --formatReadName         Default: ".fastq.gz"

    Results directory
    --resultsDir            Path to results directory, default: "./Results/"

    Assembly directory
    --assemblyDir           Path to assembly directory, default: ".Results/Assembly/"

    Script
    --renameHead            Path to rename_fasta_header.py, default: "./Tools/rename_fasta_header.py"
    --selectGetAsm          Path to script_selection_assembly.sh, default: "./Tools/script_selection_assembly.sh"
    --multi2one             Path to script convert_multiline_oneline.sh, default: "./Tools/convert_multiline_oneline.sh"

    GetOrganelle
    --getIndex             Index of GetOrganelle, default: "embplant_mt,embplant_pt"
    --getKmer              Size of kmers, default: "21,45,65,85,105"

    Sqtk
    --seqtkSubsamp         Subsampling, default: 2000000. Set to 0 in order to disable subsampling.

    FastPlast
    --fastIndex            Index of Fast-Plast, default: "Brassicales"

    OrgAsm
    --orgasmProbes         Index of ORGanelle ASeMbler, default: "protChloroArabidopsis"

    Nucmer
    --nucmerRef            Path to Fasta reference for alignment, default: "./Tools/*.fasta"

    Mummer
    --mummerAxe            Size of X-axis (fonction of genome's size), default (plastome): "'[0:154000]'"
    --mummerFormatOut      Format of the plot, default: "png"

    Mafft
    --mafftMethod          Alignment methods, default: "auto"

    Raxml
    --raxmlModel           Model uses by RAxML, default: "GTRGAMMAI"

    Each of the previous parameters can be specified as command line options or in the config file

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

process getorganelle_index {

    conda "bioconda::getorganelle"

    input:
    val index

    output:
    val("")

    script:
    """
    get_organelle_config.py -a ${params.getIndex}
    """
}

process getorganelle {

    conda "bioconda::getorganelle"

    tag "${sampleId}"

    //publishDir "${results}/Assembly", mode: 'copy'

    errorStrategy 'ignore'

    input:
    val index
    tuple val(sampleId), path(reads)

    output :
    tuple val(sampleId), val("get"), path("${sampleId}/embplant_pt.K*.complete.graph1.1.path_sequence.fasta"), path("${sampleId}/embplant_pt.K*.complete.graph1.2.path_sequence.fasta")

    script:
    """
    get_organelle_from_reads.py \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -o ${sampleId} \
    -k ${params.getKmer} \
    -F ${params.getIndex} \
    -t ${task.cpus}
    """
}

process seqtk {

    conda "bioconda::fusioncatcher-seqtk"

    tag "${sampleId}"
    
    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), path("${sampleId}_sub1.fq.gz"), path("${sampleId}_sub2.fq.gz")

    script:
    """
    seqtk sample -s666 ${reads[0]} ${params.seqtkSubsamp} |gzip -c - > ${sampleId}_sub1.fq.gz
    seqtk sample -s666 ${reads[1]} ${params.seqtkSubsamp} |gzip -c - > ${sampleId}_sub2.fq.gz
    """
}

process fastplast {

    publishDir "${results}/Assembly", mode: 'copy'

    tag "${sampleId}"

    errorStrategy 'ignore' 

    input:
    tuple val(sampleId), path(sample_R1), path(sample_R2)

    output:
    tuple val(sampleId), val("fpt"), path("${sampleId}_fpt.fasta"), emit: nucmer
    path("${sampleId}_fpt.fasta"), emit: mafft

    script:
    """
    fast-plast.pl \
    -1 ${sample_R1} \
    -2 ${sample_R2} \
    --name ${sampleId} \
    --bowtie_index ${params.fastIndex} \
    --threads ${task.cpus}

    mv "${sampleId}/Final_Assembly/${sampleId}_FULLCP.fsa" "${sampleId}_fpt.fasta"
    """
}

process orgasm {

    publishDir "${results}/Assembly", mode: 'copy'

    errorStrategy 'ignore'

    tag "${sampleId}"

    input:
    tuple val(sampleId), path(sample_R1), path(sample_R2)

    output:
    tuple val(sampleId), val("org"), path("${sampleId}_org.fasta")

    script:
    """
    oa index --estimate-length=0.9 ${sampleId} ${sample_R1} ${sample_R2}
    oa buildgraph --probes ${params.orgasmProbes} ${sampleId} ${sampleId}
    oa unfold ${sampleId} ${sampleId} > ${sampleId}_orgasm.fasta
    ${params.multi2one} "${sampleId}_orgasm.fasta" "${sampleId}_oneLine_org.fasta"
    python ${params.renameHead} -i "${sampleId}_oneLine_org.fasta" -n "${sampleId}" -o "${sampleId}_org.fasta"
    """
}

process select_assembly {

    tag "${sampleId}"

    publishDir "${results}/Assembly", mode: 'copy'

    input:
    tuple val(sampleId), val(assembler), path(assembly_1), path(assembly_2)
    
    output:
    path("${sampleId}_getGood.fasta"), emit: mafft
    tuple val(sampleId), val(assembler), path("${sampleId}_get.fasta"), emit: nucmer
    
    script:
    """
    python ${params.renameHead} -i ${assembly_1} -n "${sampleId}_get_1" -o "${sampleId}_get_1.fasta"
    python ${params.renameHead} -i ${assembly_2} -n "${sampleId}_get_2" -o "${sampleId}_get_2.fasta"
    cat "${sampleId}_get_1.fasta" "${sampleId}_get_2.fasta" > ${sampleId}_get.fasta
    bash ${params.selectGetAsm} "${sampleId}_get_1.fasta" "${sampleId}_get_2.fasta" "${sampleId}_getGood.fasta"
    """
}

process nucmer {

    conda "bioconda::mummer"

    tag "${sampleId}"

    input:
    tuple val(sampleId), val(assembler), path(assembly)

    output:
    tuple val(sampleId), val(assembler), path("${sampleId}.delta")

    script:
    """
    nucmer -p ${sampleId} ${params.nucmerRef} ${assembly}
    """
}

process mummer {

    conda "bioconda::mummer conda-forge::gnuplot"

    tag "${sampleId}"

    publishDir "${results}/Mummer", mode: 'move'

    input:
    tuple val(sampleId), val(assembler), path(align_from_nucmer)
    output:
    path("${sampleId}_${assembler}.png")

    script:
    """
    mummerplot \
    -x ${params.mummerAxe} \
    -p ${sampleId}_${assembler} \
    -t ${params.mummerFormatOut} \
    ${align_from_nucmer}
    """
}

process mafft {

    conda "bioconda::mafft"

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

process raxml {

    conda "bioconda::raxml"

    publishDir "${results}/RAxML/", mode: 'move'

    input:
    path(multi_fasta_align)

    output:
    path("RAxML_bootstrap.model_${params.raxmlModel}")

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

///////////////////////////////////////////////////////////
//////////////////////     Workflow     ///////////////////
///////////////////////////////////////////////////////////

workflow getorganelle_wf {
    take:
    data

    main:
    getorganelle_index(params.getIndex)
    getorganelle(getorganelle_index.out, data)
    select_assembly(getorganelle.out)
    nucmer(select_assembly.out.nucmer)
    mummer(nucmer.out)
}

workflow fastplast_wf {
    take:
    subsampling

    main:
    fastplast(subsampling)
    nucmer(fastplast.out.nucmer)
    mummer(nucmer.out)
}

workflow orgasm_wf {
    take:
    subsampling

    main:
    orgasm(subsampling)
    nucmer(orgasm.out)
    mummer(nucmer.out)
}

workflow all_assembler_wf {
    take:
    data

    main:
    getorganelle_wf(data)
    seqtk(data)
    fastplast_wf(seqtk.out)
    orgasm_wf(seqtk.out)
}

workflow analysing_get_wf {
    take:
    data

    main:
    getorganelle_index(params.getIndex)
    getorganelle(getorganelle_index.out, data)
    select_assembly(getorganelle.out)
    nucmer(select_assembly.out.nucmer)
    mummer(nucmer.out)
    mafft(select_assembly.out.mafft.collectFile(name: 'multi_fasta', newLine: true))
    raxml(mafft.out)
}

workflow analysing_fast_wf {
    take:
    data

    main:
    seqtk(data)
    fastplast(seqtk.out)
    nucmer(fastplast.out.nucmer)
    mummer(nucmer.out)
    mafft(fastplast.out.mafft.collectFile(name: 'multi_fasta', newLine: true))
    raxml(mafft.out)
}

workflow fromAsm_wf {
    take:
    data

    main:
    mafft(data.collectFile(name: 'multi.fasta', newLine: true))
    raxml(mafft.out)
}

///////////////////////////////////////////////////////////
//////////////////    Main Workflow     ///////////////////
///////////////////////////////////////////////////////////

workflow {

    if (params.workflow=="assembling" && params.assembler=="") { 
        data = Channel.fromFilePairs("${params.readFiles}", checkIfExists:true)
        all_assembler_wf(data) }

    else if (params.workflow=="assembling" && params.assembler=='getorganelle') { 
        data = Channel.fromFilePairs("${params.readFiles}", checkIfExists:true)
        getorganelle_wf(data) }

    else if (params.workflow=="assembling" && params.assembler=='fastplast') { 
        data = Channel.fromFilePairs("${params.readFiles}", checkIfExists:true)
        fastplast_wf(data) }

    else if (params.workflow=="assembling" && params.assembler=='orgasm') { 
        data = Channel.fromFilePairs("${params.readFiles}", checkIfExists:true)
        orgasm_wf(data) }  

    else if (params.workflow=="analyzing" && (params.assembler=="getorganelle" || params.assembler=="")) { 
        data = Channel.fromFilePairs("${params.readFiles}", checkIfExists:true)
        analysing_get_wf(data) }

    else if (params.workflow=="analyzing" && params.assembler=="fastplast") { 
        data = Channel.fromFilePairs("${params.readFiles}", checkIfExists:true)
        analysing_fast_wf(data) }

    else if (params.workflow=="fromAsm") { 
        data = Channel.fromPath("${params.assemblyDir}")
        fromAsm_wf(data) }

    else
    { println """
    
    Error : Any workflow selected or unknown assembler
    
    """
    helpMessage()
    exit 0 }
}

workflow.onComplete{println("Workflow execution completed sucessfully!")}