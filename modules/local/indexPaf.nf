process indexPaf{
    input:
        file ref
    output:
        path 'ref.mmi', emit: ref_idx
        path '*.fa', emit: ref_unzipped
    script:
        """
        minimap2 -d ref.mmi $ref
        if [[ $ref == *.gz ]]
        then
            gunzip $ref -c > ${ref.getSimpleName()}.fa
        else
            cat $ref > ${ref.getSimpleName()}_unzipped.fa
        fi
        """
}