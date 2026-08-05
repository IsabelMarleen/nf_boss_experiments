process getMapRef {
    input:
        val ref_link
    output:
        path '*.fa', emit: ref_path
    script:
    """
    wget -O ref_for_mapping.fa ${ref_link}
    """

}