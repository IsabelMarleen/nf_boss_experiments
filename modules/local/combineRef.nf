process combineRef{
    input:
        path refs
    output:
        path '*.fa.gz', emit: ref_path
    script:
    ref_base = "${refs.first().getBaseName(3)}"
    """
    zcat ${refs} | awk '{gsub(/^>/,">chr");}1' | gzip > "${ref_base}.${params.chromosomes.join(".")}.fa.gz"
    """
}