process getFile {
    input:
        val link
    output:
        path '*.{fa.gz,bed,tar.gz}', emit: path
    script:
    """
    wget "${link}"
    """

}