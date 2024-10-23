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