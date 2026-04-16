// Map paf trunc Step 4 of prepare input for BRsim
process mapPaf_trunc {
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1'
    time {30.min * params.readnumber.values().size() * params.chromosomes.size() * task.attempt}
    publishDir "${params.output_dir_br_input}", mode: 'symlink'
    conda "${params.conda_base_dir}/boss2"
    input: 
        path aligned_fastq_trunc
        path ref
    output: 
        path "*.paf", emit: mappings_trunc

    script:
    """
    # from shutil import which

    # mm2 = which("minimap2", path='/'.join(executable.split('/')[0:-1])).strip()
    minimap2 -x map-ont -t 32 --secondary=no -c $ref $aligned_fastq_trunc > ${aligned_fastq_trunc.getSimpleName()}.paf
    """
}