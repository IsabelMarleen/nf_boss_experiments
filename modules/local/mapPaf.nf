process mapPaf {
    tag "${input.getSimpleName()}"
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    conda "${params.conda_envs}/simulation.yaml"
    input:
        file input
        file index
    output:
        path '*.paf', emit: ref_path
    script:
    """
    minimap2 -x map-ont -t 32 --secondary=no -c $index $input -o ${input.getSimpleName()}.paf
    """

}