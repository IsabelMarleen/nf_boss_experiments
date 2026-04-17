// Scan fq offsets Step 2 of prepare input for BRsim
process scan_offsetsFq {
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    publishDir "${params.output_dir_br_input}", mode: 'symlink'
    conda "${params.conda_base_dir}/boss2"
    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/99/99949ae2aa4e34667e718a42be92c9064fdd3dae10afeb04fd99b3162c3adfa5/data"
    input: 
        path aligned_fastq
        path sampler_func_file, stageAs: 'sampler.py'

    output: 
        path "*.offsets.npy", emit: fq_offsets

    script:
    """
    scan_offsets_fq_nf.py $aligned_fastq
    """
}