process getRef {
    clusterOptions '--mem=2G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_vi}", mode: 'symlink'
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