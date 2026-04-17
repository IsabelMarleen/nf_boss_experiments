process combineRef{
    clusterOptions '--mem=2G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_vi}", mode: 'symlink'
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