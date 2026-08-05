process combineBarcodedFq{
    input: 
        path barcoded_fastq
    output: 
        path "*_barcoded.fq", emit:combined_barcoded_fq

    script:
    """
    cat ${barcoded_fastq} | \
    awk '{OFS="\\t"; getline seq; \
                getline sep; \
                getline qual; \
                print \$0,seq,sep,qual}' | \
    shuf | \
    awk -F"\\t" '{OFS="\\n"; print \$1,\$2,\$3,\$4}' \
    > ${params.exp_name}_combined_barcoded.fq
    """
}