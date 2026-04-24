process pileup {
    tag "${meta.join("_")}"
    input:
        tuple val(meta), path(bam)
        tuple val(meta2), path(bai)
        file otu_file
    output: 
        path 'tmp.csv', emit: pup
        path "${meta.join("_")}.csv", emit: csv
    script:
    if (params.barcodes == null){
        otu = meta.get(1)
        assert otu == meta2.get(1)
    }
    else{
        otu = meta.get(2)
        assert otu == meta2.get(2)
    }
    """
    samtools mpileup -Q 0 ${bam} > ${meta.join("_")}.pup &&
    process_pileup.py ${meta.join("_")}.pup $otu ${meta.join("_")}.csv ${bam.size()}
    """
}