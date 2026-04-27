process addBarcodes{
    input: 
        path simulated_fastq
    output: 
        path "*_barcoded.fq", emit:barcoded_fq

    script:
    labels = params.barcodes.keySet().join("|")
    def matcher = simulated_fastq =~ "$labels"
    barcode = params.barcodes[matcher[0]]
    """
    sed '1~4 s/\$/ barcode=$barcode/' ${simulated_fastq} > ${simulated_fastq.getSimpleName()}_barcoded.fq
    """
}