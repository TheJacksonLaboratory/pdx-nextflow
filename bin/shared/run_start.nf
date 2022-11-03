process RUN_START {
    publishDir "${params.pubdir}", mode:'copy'

    output:
    file("pipeline_running.txt")

    script:
    """
    touch pipeline_running.txt
    """
}