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

    Command : nextflow run ChloroBras.nf --workflow [assembling/analyzing/fromAsm]

    REQUIRED parameter

    -profile [standard]/slurm,      Select profile standard (local) or slurm. Default: standard          
             singularity/conda      Select profile singularity or conda. (FastPlast and Orgasm are only available with singularity, even in conda profile)
                                                                         (Mummer is only available with conda, even in singularity profile)

    --workflow [assembling/analyzing/fromAsm]       assembling : assembles genomes with the three assemblers, allows quality assessment via dotplot
                                                    analyzing : assemble genomes with [getorganelle] or Fastplast and create phylogenetic tree
                                                    fromAsm : mafft alignement and Raxml tree from assemblies in ./Results/Assembly/
    
    OPTIONAL parameter

    Assembler
    --assembler             Choose assembler to use (getorganelle, fastplast or orgasm), default: all for assembling workflow
                                                                                         default: getorganelle for analysing workflow
    
    Reads directory
    --readDir                Default: "./Data"
    --baseReadName           Default: "_R{1,2}"     ex: name_R1.fastq.gz & name_R2.fastq.gz
    --formatReadName         Default: ".fastq.gz"

    Results directory
    --resultsDir            Path to results directory, default: "./Results/"

    Assembly directory
    --assemblyDir           Path to assembly directory, default: "./Results/Assembly/"
    --formatAsm             Default: ".fasta"

    GetOrganelle
    --getIndex             Index of GetOrganelle, default: "embplant_mt,embplant_pt"
    --getKmer              Size of kmers, default: "21,45,65,85,105"

    Sqtk
    --seqtkSubsamp         Subsampling, default: 2000000.

    FastPlast
    --fastIndex            Index of Fast-Plast, default: "Brassicales"

    OrgAsm
    --orgasmProbes         Index of ORGanelle ASeMbler, default: "protChloroArabidopsis"

    Mummer
    --nucmerRef            Path to Fasta reference for alignment, default: "./Data/brassica_oleracea.fasta"
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

    label 'getorganelle'

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

    label 'getorganelle' 
    label 'process_high'

    tag "${sampleId}"

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

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
    tuple val(sampleId), path("${sampleId}_sub1.fq.gz"), path("${sampleId}_sub2.fq.gz")

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

    input:
    tuple val(sampleId), path(sample_R1), path(sample_R2)

    output:
    tuple val(sampleId), val("fpt"), path("${sampleId}_fpt.fasta"), emit: mummer
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

    label 'process_high'

    publishDir "${results}/Assembly/", mode: 'copy'

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
    convert_multiline_oneline.sh "${sampleId}_orgasm.fasta" "${sampleId}_oneLine_org.fasta"
    rename_fasta_header.py -i "${sampleId}_oneLine_org.fasta" -n "${sampleId}" -o "${sampleId}_org.fasta"
    """
}

process mummer {

    tag "${sampleId}"

    publishDir "${results}/Mummer", mode: 'move'

    input:
    tuple val(sampleId), val(assembler), path(assembly)

    output:
    path("${sampleId}_${assembler}.png")

    script:
    """
    nucmer -p ${sampleId} ${params.nucmerRef} ${assembly}

    mummerplot \
    -x ${params.mummerAxe} \
    -p ${sampleId}_${assembler} \
    -t ${params.mummerFormatOut} \
    "${sampleId}.delta"
    """
}

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

process raxml {

    label 'process_high'

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
//////////////////////     Sub-Workflow     ///////////////
///////////////////////////////////////////////////////////

workflow getorganelle_wf {
    take:
    data

    main:
    getorganelle_index(params.getIndex)
    getorganelle(getorganelle_index.out, data)
    mummer(getorganelle.out.mummer)
}

workflow fastplast_wf {
    take:
    subsampling

    main:
    fastplast(subsampling)
    mummer(fastplast.out.mummer)
}

workflow orgasm_wf {
    take:
    subsampling

    main:
    orgasm(subsampling)
    mummer(orgasm.out)
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
    mummer(getorganelle.out.mummer)
    mafft(getorganelle.out.mafft.collectFile(name: 'multi_fasta', newLine: true))
    raxml(mafft.out)
}

workflow analysing_fast_wf {
    take:
    data

    main:
    seqtk(data)
    fastplast(seqtk.out)
    mummer(fastplast.out.mummer)
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

///////////////////////////////////////////////////////
//////////////////     Workflow     ///////////////////
///////////////////////////////////////////////////////

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
        data = Channel.fromPath("${params.asmFiles}")
        fromAsm_wf(data) }

    else
    { println """
    
    Error : Any workflow selected or unknown assembler
    
    """
    helpMessage()
    exit 0 }
}

workflow.onComplete{println("Workflow execution completed sucessfully!")}