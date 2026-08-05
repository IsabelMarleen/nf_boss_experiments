process pileup {
    tag "${meta.id}"
    input:
        tuple val(meta), path(bam), path(bai)
        file otu_file
    output: 
        path "${meta.id}.csv", emit: csv
    script:
    otu = bam.simpleName.split(/_/)[-2]
    """
    samtools mpileup -Q 0 ${bam} > ${meta.id}.pup &&
    process_pileup.py ${meta.id}.pup $otu ${meta.id}.csv
    """
}