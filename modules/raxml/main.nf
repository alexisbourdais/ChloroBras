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