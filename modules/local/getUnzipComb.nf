process getUnzipComb {
    input:
        val link
    output:
        path '*.fa', emit: ref_path
    script:
    """
    wget -O folder.zip ${link}
    unzip folder.zip
    cat */Genomes/* > reference.fa
    """

}