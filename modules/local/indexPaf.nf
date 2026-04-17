process indexPaf{
    clusterOptions '--mem=8G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:10:00'
    conda "${params.conda_envs}/simulation.yaml"
    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f666f3a7c3db2bc3beead920385828167c735d1b3e47530303c12c0cd81329a0/data"
    input:
        file ref
    output:
        path 'ref.mmi', emit: ref_idx
        path '*.fa', emit: ref_unzipped
    script:
        """
        minimap2 -d ref.mmi $ref
        gunzip $ref -c > ${ref.getSimpleName()}.fa
        """
}