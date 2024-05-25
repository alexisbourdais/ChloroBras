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

    Each of the following parameters can be specified as command line options or in the config file

    Reads directory
    --readsFiles                    Path to input data, default: "${params.baseDir}/Samples/*_R{1,2}.fastq.gz"

    Results directory
    --resultsDir                    Path to results directory, default: "${params.baseDir}/Results/"

    GetOrganelle
    --getorganelle_index            Index of GetOrganelle, default: "embplant_pt"
    --getorganelle_kmer             Size of kmers, default: "21,45,65,85,105"

    Rename_headers
    --rename_script                 Path to rename_fasta_header.py, default: "${params.baseDir}/Tools/rename_fasta_header.py"

    Nucmer
    --nucmer_ref                    Path to Fasta reference for alignment, default: "${params.baseDir}/Tools/brassica_oleracea.fasta"

    Mummer
    --mummer_axe                    Size of X-axis (fonction of genome's size), default (plastome): "'[0:154000]'"
    --mummer_format_output          Format of the plot, default: "png"

    Select_assembly
    --select_assembly_script        Path to script_selection_assembly.sh, default: "${params.baseDir}/Tools/script_selection_assembly.sh"

    Mafft
    --mafft_method                  Alignment methods, default: "auto"

    Raxml
    --raxml_model                   Model uses by RAxML, default: "GTRGAMMAI"

    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
        helpMessage()
        exit 0
}

///////////////////////////////////////////////////////////
//////////////////////     Processes    ///////////////////
///////////////////////////////////////////////////////////

results = file(params.resultsDir)

process getorganelle_index {

    conda "${params.baseDir}/Tools/getorganelle_env.yml"

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

    conda "${params.baseDir}/Tools/getorganelle_env.yml"

    tag "${sampleId}"

    input:
    val index
    tuple val(sampleId), path(reads)

    output :
    tuple val(sampleId), path("${sampleId}/embplant_pt.K*.complete.graph1.1.path_sequence.fasta"), path("${sampleId}/embplant_pt.K*.complete.graph1.2.path_sequence.fasta")

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

process rename_headers {

    tag "${sampleId}"

    input:
    tuple val(sampleId), path(assembly_getorganelle_1), path(assembly_getorganelle_2)
    
    output:
    tuple val(sampleId), path("${sampleId}_get_1"), path("${sampleId}_get_2")
    
    script:
    """
    python ${params.rename_script} -i ${assembly_getorganelle_1} -n "${sampleId}_get_1" -o "${sampleId}_get_1"
    python ${params.rename_script} -i ${assembly_getorganelle_2} -n "${sampleId}_get_2" -o "${sampleId}_get_2"
    """
}

process concaten {

    tag "${sampleId}"

    input:
    tuple val(sampleId), path(assembly_getorganelle_1), path(assembly_getorganelle_2)

    output:
    tuple val(sampleId), path("${sampleId}.fasta")

    script:
    """
    cat ${assembly_getorganelle_1} ${assembly_getorganelle_2} > ${sampleId}.fasta
    """
}

process nucmer {

    conda "${params.baseDir}/Tools/mummer_env.yml"

    tag "${sampleId}"

    input:
    tuple val(sampleId), path(assembly_cat)

    output:
    tuple val(sampleId), path("${sampleId}.delta")

    script:
    """
    nucmer -p ${sampleId} ${params.nucmer_ref} ${assembly_cat}
    """
}

process mummer {

    conda "${params.baseDir}/Tools/mummer_env.yml"

    tag "${sampleId}"

    publishDir "${results}/${sampleId}/", mode: 'move'

    input:
    tuple val(sampleId), path(align_from_nucmer)

    output:
    tuple val(sampleId), path("${sampleId}_get.png")

    script:
    """
    mummerplot \
    -x ${params.mummer_axe} \
    -p ${sampleId}_get \
    -t ${params.mummer_format_output} \
    ${align_from_nucmer}
    """
}

process select_assembly {

    tag "${sampleId}"

    input:
    tuple val(sampleId), path(assembly_getorganelle_1), path(assembly_getorganelle_2)

    output:
    path("${sampleId}_graph_right")

    script:
    """
    bash ${params.select_assembly_script} ${assembly_getorganelle_1} ${assembly_getorganelle_2} ${sampleId}_graph_right
    """
}

process mafft {

    conda "${params.baseDir}/Tools/mafft_env.yml"

    tag "Alignment"

    publishDir "${results}/MAFFT", mode: 'copy'

    input:
    path(multi_fasta)

    output:
    path("multi_fasta_align.fasta")

    script:
    """
    mafft --${params.mafft_method} ${multi_fasta} > multi_fasta_align.fasta
    sed -i -e 's/_get_1_1//g' -e 's/_get_2_1//g' multi_fasta_align.fasta
    """
}

process raxml {

    conda "${params.baseDir}/Tools/raxml_env.yml"

    tag "Creation of the phylogenetic tree"

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

workflow {
    data = Channel.fromFilePairs("${params.readsFiles}", checkIfExists:true)
    getorganelle_index(params.getorganelle_index)
    getorganelle(getorganelle_index.out, data)
    rename_headers(getorganelle.out)
    concaten(rename_headers.out)
    nucmer(concaten.out)
    mummer(nucmer.out)
    select_assembly(rename_headers.out)
    mafft(select_assembly.out.collectFile(name: 'multi_fasta', newLine: true))
    raxml(mafft.out)
}

workflow.onComplete{println("Workflow execution completed sucessfully!")}