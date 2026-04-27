process separateTarget {
    tag "${paf.getSimpleName()}"
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