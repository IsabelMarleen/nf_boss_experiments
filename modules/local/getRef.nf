process getRef {
    input:
        val ref_link
        val chr
    output:
        path '*.fa.gz', emit: ref_path
    script:
    """
    wget "${ref_link}${chr}.fa.gz"
    """

}