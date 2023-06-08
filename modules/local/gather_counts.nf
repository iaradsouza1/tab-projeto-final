process GATHER_COUNTS {
    tag "$meta.id"
    label "process_medium"

    container "rocker/tidyverse:4.2.2"

    input:
    path ("subread/*.featureCounts.txt")

    output:
    path "count_table.txt"
    path "count_table.rds"

    when:
    task.ext.when == null || task.ext.when

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in iaradsouza1/tab-projeto-final/bin
    """
    get_counts.R \\
        subread

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """


}
