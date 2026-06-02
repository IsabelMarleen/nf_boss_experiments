process combineLargeBarcodedFq{
    input: 
        path barcoded_fastq
    output: 
        path "*_barcoded.fq", emit:combined_barcoded_fq

    script:
    if (params.seed == null){
        seed = 27
    }
    else{
        seed = params.seed
    }
    """
    combine_barcoded.py ${params.exp_name}_combined_barcoded.fq ${params.readnumber.values().sum()*0.6} ${barcoded_fastq} ${seed}
    """
}