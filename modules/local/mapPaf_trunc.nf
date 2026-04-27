// Map paf trunc Step 4 of prepare input for BRsim
process mapPaf_trunc {
    input: 
        path aligned_fastq_trunc
        path ref
    output: 
        path "*.paf", emit: mappings_trunc

    script:
    """
    minimap2 -x map-ont -t 32 --secondary=no -c $ref $aligned_fastq_trunc > ${aligned_fastq_trunc.getSimpleName()}.paf
    """
}