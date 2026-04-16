// Scan paf offsets Step 5 of prepare input for BRsim
process scan_offsets_Paf {
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    publishDir "${params.output_dir_br_input}", mode: 'symlink'
    conda "${params.conda_base_dir}/boss2"
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