process getConsensus {
    input:
        tuple path(vcf), path(vcf_idx)
        path ref
    output:
        path '*.fa', emit:consensus
    script:
    """
    bcftools consensus -f $ref $vcf -o consensus_${vcf.getBaseName(2)}.fa
    """
}