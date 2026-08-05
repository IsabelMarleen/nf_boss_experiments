process addBarcodes{
    input: 
        path simulated_fastq
        val barcoding_needed
    output: 
        path "*_barcoded.fq", emit:barcoded_fq
    when:
        barcoding_needed
    script:
    labels = params.barcodes.keySet().join("|")
    def matcher = simulated_fastq =~ "$labels"
    barcode = params.barcodes[matcher[0]]
    """
    sed '1~4 s/\$/ barcode=$barcode/' ${simulated_fastq} > ${simulated_fastq.getSimpleName()}_barcoded.fq
    """
}