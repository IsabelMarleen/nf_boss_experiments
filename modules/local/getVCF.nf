process getVCF {
    clusterOptions '--mem=2G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_vi}", mode: 'symlink'
    input:
        tuple val(vcf_link), val(idx_link)
    output:
        tuple path('*.vcf.gz'), path('*.vcf.gz.tbi')
    script:
    """
    wget $vcf_link
    wget $idx_link
    """

}