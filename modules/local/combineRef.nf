process combineRef{
    input:
        path refs
    output:
        path '*.fa.gz', emit: ref_path
    script:
    ref_base = "${refs.first().getBaseName(3)}"
    """
    cat ${refs} > "${ref_base}.${params.chromosomes.join(".")}.fa.gz"
    """
}