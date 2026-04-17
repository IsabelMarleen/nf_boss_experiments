process separateTarget {
    tag "${paf.getSimpleName()}"
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:10:00'
    conda "${params.conda_envs}/simulation.yaml"
    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/33/33d3e3989ad86b8390397b202c6db499cae435dc34174d67fb36027adf25f41c/data"
    input:
        tuple file(paf), file(reads)
        file otu
    output:
        path '*.fq', emit: fq
    script:
    """
    separate_by_target.py $paf $reads
    """

}