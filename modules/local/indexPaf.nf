process indexPaf{
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