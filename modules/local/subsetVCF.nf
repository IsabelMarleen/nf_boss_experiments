process subsetVCF {
    input:
        tuple path(full_vcf), path(full_idx)
    output:
        tuple path('*.vcf.gz'), path('*.vcf.gz.tbi')

    script:
    if (params.chromosomes.size() > 1){
        chr_string = params.chromosomes.collect{"chr" + it}.join(",")
        chr_prefix = params.chromosomes.join("_")
    }
    else{
        chr_string = "chr${params.chromosomes.first()}"
        chr_prefix = params.chromosomes.first()
    }
    """
    bcftools view -r "$chr_string" -Oz $full_vcf > "${chr_prefix}_${full_vcf}"
    bcftools index -t "${chr_prefix}_${full_vcf}"
    """
}