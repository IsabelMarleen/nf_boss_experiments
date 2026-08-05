process extractTar {
    input:
        path directory
        val strat_file
    output:
        path "*[!.tsv]", emit: path
        path "${strat_file}", emit: strat
    script:
    """
    tar -xf $directory
    touch $strat_file
    """

}