process separateTarget {
    tag "${paf.getSimpleName()}"
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:10:00'
    conda "${params.conda_envs}/simulation.yaml"
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