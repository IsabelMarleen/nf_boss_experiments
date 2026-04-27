// Scan paf offsets Step 5 of prepare input for BRsim
process scan_offsets_Paf {
    input: 
        path mappings
        path mappings_trunc
    output: 
        path "${mappings.getSimpleName()}.paf.offsets", emit: mappings_offsets
        path "*_trunc.paf.offsets", emit: mappings_offsets_trunc

    script:
    """
    scan_offsets_paf_nf.py ${mappings} ${mappings_trunc}
    """
}