process recordUnblocks {
    tag "${fq.getSimpleName()}"
    clusterOptions '--mem=24G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:10:00'
    conda "${params.conda_envs}/simulation.yaml"
    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/82/82850b60e32398db58f09fd538b44ac29fd7b64f62adc1661477614dfbd85066/data"
    input:
        path fq
    output:
        path '*.csv', emit: csv
    script:
    """
    record_unblocks.py $fq ${params.mu} > ${fq.getSimpleName()}_unblocks.csv
    """
}