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


    Command : nextflow run ChloroBras.nf --workflow [test/analysis]


    REQUIRED parameter

    Workflow
    --workflow [test/analysis]      test : assembles genomes with the three assemblers, allows quality assessment via dotplot
                                    analysis : assemble genomes with getorganelle and create phylogenetic tree

    OPTIONAL parameter

    Executor
    --executor                      Choose the executor (local or slurm). Default: local

    Analysis assembler
    --analysis_assembler [getorganelle/fastplast] Change the assembler used to the analysis workflow. Default: getorganelle

    Reads directory
    --readsFiles                    Path to input data, default: "${params.baseDir}/Samples/*_R{1,2}.fastq.gz"

    Results directory
    --resultsDir                    Path to results directory, default: "${params.baseDir}/Results/"

    GetOrganelle
    --getorganelle_index            Index of GetOrganelle, default: "embplant_pt"
    --getorganelle_kmer             Size of kmers, default: "21,45,65,85,105"

    Sqtk
    --seqtk_nb_read                 Number of reads to keep, default: 2000000

    FastPlast
    --fastplast_index               Index of Fast-Plast, default: "Brassicales"

    OrgAsm
    --orgasm_probes                 Index of ORGanelle ASeMbler, default: "protChloroArabidopsis"

    Rename_headers

    Nucmer
    --nucmer_ref                    Path to Fasta reference for alignment, default: "${params.baseDir}/Tools/brassica_oleracea.fasta"

    Mummer
    --mummer_axe                    Size of X-axis (fonction of genome's size), default (plastome): "'[0:154000]'"
    --mummer_format_output          Format of the plot, default: "png"

    Select_assembly

    Mafft
    --mafft_method                  Alignment methods, default: "auto"

    Raxml
    --raxml_model                   Model uses by RAxML, default: "GTRGAMMAI"

    Script
    --rename_script                 Path to rename_fasta_header.py, default: "${params.baseDir}/Tools/rename_fasta_header.py"
    --select_assembly_script        Path to script_selection_assembly.sh, default: "${params.baseDir}/Tools/script_selection_assembly.sh"
    --multi2one_script              Path to script convert_multiline_oneline.sh, default: "${params.baseDir}/Tools/convert_multiline_oneline.sh


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
    get_organelle_config.py --add ${index}
    """
}

process getorganelle {

    conda "bioconda::getorganelle"

    tag "${sampleId}"

    //publishDir "${results}/${sampleId}/GetOrganelle", mode: 'copy'

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
    -k ${params.getorganelle_kmer} \
    -F ${params.getorganelle_index} \
    -t 1
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
    seqtk sample -s666 ${reads[0]} ${params.seqtk_nb_read} |gzip -c - > ${sampleId}_sub1.fq.gz
    seqtk sample -s666 ${reads[1]} ${params.seqtk_nb_read} |gzip -c - > ${sampleId}_sub2.fq.gz
    """
}

process fastplast {

    publishDir "${results}/${sampleId}/Fastplast", mode: 'copy'

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
    --bowtie_index ${params.fastplast_index} \

    mv "${sampleId}/Final_Assembly/${sampleId}_FULLCP.fsa" "${sampleId}_fpt.fasta"
    """
}

process orgasm_index {

    label 'orgasm'

    tag "${sampleId}"

    errorStrategy 'ignore'

    input:
    tuple val(sampleId), path(sample_R1), path(sample_R2)

    output:
    tuple val(sampleId), path("${sampleId}.odx")

    script:
    """
    oa index --estimate-length=0.9 ${sampleId} ${sample_R1} ${sample_R2}
    """
}

process orgasm_buildgraph {

    label 'orgasm'

    tag "${sampleId}"

    errorStrategy 'ignore'

    input:
    tuple val(sampleId), path(index)

    output:
    tuple val(sampleId), path(index), path("${sampleId}.oas")

    script:
    """
    oa buildgraph --probes ${params.orgasm_probes} ${sampleId} ${sampleId}
    """
}

process orgasm_unfold {

    label 'orgasm'

    tag "${sampleId}"

    publishDir "${results}/${sampleId}/Orgasm/", mode: 'copy'

    errorStrategy 'ignore'

    input:
    tuple val(sampleId), path(index), path(graph)

    output:
    tuple val(sampleId), val("org"), path("${sampleId}_org.fasta")

    script:
    """
    oa unfold ${sampleId} ${sampleId} > ${sampleId}_orgasm.fasta
    ${params.multi2one} "${sampleId}_orgasm.fasta" "${sampleId}_oneLine_org.fasta"
    python ${params.rename_script} -i "${sampleId}_oneLine_org.fasta" -n "${sampleId}" -o "${sampleId}_org.fasta"
    """
}

process select_assembly {

    tag "${sampleId}"

    publishDir "${results}/${sampleId}/GetOrganelle/", mode: 'copy'

    input:
    tuple val(sampleId), val(assembler), path(assembly_1), path(assembly_2)
    
    output:
    path("${sampleId}_getGood.fasta"), emit: mafft
    tuple val(sampleId), val(assembler), path("${sampleId}_get.fasta"), emit: nucmer
    
    script:
    """
    python ${params.rename_script} -i ${assembly_1} -n "${sampleId}_get_1" -o "${sampleId}_get_1.fasta"
    python ${params.rename_script} -i ${assembly_2} -n "${sampleId}_get_2" -o "${sampleId}_get_2.fasta"
    cat "${sampleId}_get_1.fasta" "${sampleId}_get_2.fasta" > ${sampleId}_get.fasta
    bash ${params.select_assembly_script} "${sampleId}_get_1.fasta" "${sampleId}_get_2.fasta" "${sampleId}_getGood.fasta"
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
    nucmer -p ${sampleId} ${params.nucmer_ref} ${assembly}
    """
}

process mummer {

    conda "bioconda::mummer conda-forge::gnuplot"

    tag "${sampleId}"

    publishDir "${results}/${sampleId}/Mummer", mode: 'move'

    input:
    tuple val(sampleId), val(assembler), path(align_from_nucmer)
    output:
    path("${sampleId}_${assembler}.png")

    script:
    """
    mummerplot \
    -x ${params.mummer_axe} \
    -p ${sampleId}_${assembler} \
    -t ${params.mummer_format_output} \
    ${align_from_nucmer}
    """
}

process mafft {

    conda "bioconda::mafft"

    tag "Alignment"

    publishDir "${results}/Mafft-Alignment/", mode: 'copy'

    input:
    path(multi_fasta)

    output:
    path("multi_fasta_align.fasta")

    script:
    """
    mafft --${params.mafft_method} ${multi_fasta} > multi_fasta_align.fasta
    sed -i -e 's/_get_1//g' -e 's/_get_2//g' multi_fasta_align.fasta
    """
}

process raxml {

    conda "bioconda::raxml"

    tag "Phylogenetic tree"

    publishDir "${results}/RAxML/", mode: 'move'

    input:
    path(multi_fasta_align)

    output:
    path("RAxML_bootstrap.model_${params.raxml_model}")

    script:
    """
    raxmlHPC \
    -s ${multi_fasta_align} \
    -n model_${params.raxml_model} \
    -d \
    -m ${params.raxml_model} \
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
    getorganelle_index(params.getorganelle_index)
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
    orgasm_index(subsampling)
    orgasm_buildgraph(orgasm_index.out)
    orgasm_unfold(orgasm_buildgraph.out)
    nucmer(orgasm_unfold.out)
    mummer(nucmer.out)
}

workflow test_assembler_wf {
    take:
    data

    main:
    getorganelle_wf(data)
    seqtk(data)
    fastplast_wf(seqtk.out)
    orgasm_wf(seqtk.out)
}

workflow analysis_get_wf {
    take:
    data

    main:
    getorganelle_index(params.getorganelle_index)
    getorganelle(getorganelle_index.out, data)
    select_assembly(getorganelle.out)
    nucmer(select_assembly.out.nucmer)
    mummer(nucmer.out)
    mafft(select_assembly.out.mafft.collectFile(name: 'multi_fasta', newLine: true))
    raxml(mafft.out)
}

workflow analysis_fast_wf {
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

///////////////////////////////////////////////////////////
//////////////////    Main Workflow     ///////////////////
///////////////////////////////////////////////////////////

workflow {

    data = Channel.fromFilePairs("${params.readsFiles}", checkIfExists:true)

    if (params.workflow=="test")
    { test_assembler_wf(data) }

    else if (params.workflow=="analysis" && params.analysis_assembler=="getorganelle")
    { analysis_get_wf(data) }

    else if (params.workflow=="analysis" && params.analysis_assembler=="fastplast")
    { analysis_fast_wf(data) }

    else
    { println """
    
    Error : Any workflow selected or unknown assembler
    
    """
    helpMessage()
    exit 0 }
}

workflow.onComplete{println("Workflow execution completed sucessfully!")}