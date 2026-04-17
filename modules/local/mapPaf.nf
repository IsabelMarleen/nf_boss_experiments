process mapPaf {
    tag "${input.getSimpleName()}"
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    conda "${params.conda_envs}/simulation.yaml"
    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f666f3a7c3db2bc3beead920385828167c735d1b3e47530303c12c0cd81329a0/data"
    input:
        file input
        file index
    output:
        path '*.paf', emit: mappings
    script:
    """
    minimap2 -x map-ont -t 32 --secondary=no -c $index $input -o ${input.getSimpleName()}.paf
    """

}