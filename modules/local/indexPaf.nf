process indexPaf{
    clusterOptions '--mem=8G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:10:00'
    conda "${params.conda_envs}/simulation.yaml"
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