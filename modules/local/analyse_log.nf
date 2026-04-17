// Analyse read log
process analyse_log {
    clusterOptions '--mem=8G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_debug}/results", mode: 'symlink'
    conda "${params.conda_envs}/simulation.yaml"
    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/8776235a94a745da74c053f4c2f0c7452be44babb551410832b1dea0a3bc086d/data"
    input: 
        path log_file
    output: 
        path "*.csv", emit: analysed_log

    script:
    """
    analyse_boss_log.py $log_file ${log_file.getSimpleName()}_log_processed.csv
    """
}