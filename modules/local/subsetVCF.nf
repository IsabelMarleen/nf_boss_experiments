process subsetVCF {
    input:
        tuple path(full_vcf), path(full_idx)
    output:
        tuple path('*.vcf.gz'), path('*.vcf.gz.tbi')

    script:
    if (params.chromosomes.size() > 1){
        chr_string = params.chromosomes.collect{"chr" + it}.join(",")
        chr_name_match = params.chromosomes.collect{"chr" + it + " " + it}.join("\n")
        chr_prefix = params.chromosomes.join("_")
    }
    else{
        chr_string = "chr${params.chromosomes.first()}"
        chr_name_match = "chr${params.chromosomes.first()} ${params.chromosomes.first()}"
        chr_prefix = params.chromosomes.first()
    }
    """
    echo "${chr_name_match}" > chr_name_match.txt
    bcftools view -r "$chr_string" -Oz $full_vcf | bcftools annotate -Oz --rename-chrs chr_name_match.txt -o "${chr_prefix}_${full_vcf}"
    bcftools index -t "${chr_prefix}_${full_vcf}"
    """
}