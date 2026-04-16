process recordUnblocks {
    tag "${fq.getSimpleName()}"
    clusterOptions '--mem=24G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:10:00'
    conda "${params.conda_envs}/simulation.yaml"
    input:
        path fq
    output:
        path '*.csv', emit: csv
    script:
    """
    record_unblocks.py $fq ${params.mu} > ${fq.getSimpleName()}_unblocks.csv
    """
}