// Truncate fq Step 1 of prepare input for BRsim
process truncateFq {
    input: 
        path aligned_fastq
    output: 
        path "*_trunc.fq", emit:trunc_fq
        path aligned_fastq, emit:full_fq

    script:
    """
    cut -c -${params.mu} ${aligned_fastq} > ${aligned_fastq.getSimpleName()}_trunc.fq
    """
}