process recordUnblocks {
    tag "${fq.getSimpleName()}"
 input:
        path fq
    output:
        path '*.csv', emit: csv
    script:
    """
    record_unblocks.py $fq ${params.mu} > ${fq.getSimpleName()}_unblocks.csv
    """
}