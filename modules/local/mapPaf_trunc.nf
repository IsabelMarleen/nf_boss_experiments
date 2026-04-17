// Map paf trunc Step 4 of prepare input for BRsim
process mapPaf_trunc {
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1'
    time {30.min * params.readnumber.values().size() * params.chromosomes.size() * task.attempt}
    publishDir "${params.output_dir_br_input}", mode: 'symlink'
    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f666f3a7c3db2bc3beead920385828167c735d1b3e47530303c12c0cd81329a0/data"
    conda "${moduleDir}/../../envs/simulation.yaml"
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