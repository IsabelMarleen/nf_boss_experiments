process mapPaf {
    tag "${input.getSimpleName()}"
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