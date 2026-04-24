// Scan fq offsets Step 2 of prepare input for BRsim
process scan_offsetsFq {
    input: 
        path aligned_fastq

    output: 
        path "*.offsets.npy", emit: fq_offsets

    script:
    """
    scan_offsets_fq_nf.py $aligned_fastq
    """
}