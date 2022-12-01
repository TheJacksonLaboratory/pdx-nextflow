process RUN_START {
    cpus 1
    memory 10.MB
    time '00:03:00'

    publishDir "${params.pubdir}", mode:'copy'

    output:
    file("pipeline_running.txt")

    script:
    """
    touch pipeline_running.txt
    """
}