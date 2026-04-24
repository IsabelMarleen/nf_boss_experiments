process getVCF {
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