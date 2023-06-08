process GATHER_COUNTS {
    label "process_medium"

    container "biocontainers/r-tidyverse:1.2.1"

    input:
    path feature_counts

    output:
    path "count_table.txt"
    path "count_table.rds"

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in iaradsouza1/tab-projeto-final/bin
    """
    get_counts.r \\
        $feature_counts

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """

}
