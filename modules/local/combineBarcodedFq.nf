process combineBarcodedFq{
    clusterOptions '--nodes=1 --cpus-per-task=8 --ntasks=1'
    memory { 64.GB * params.readnumber.values().size() * task.attempt }
    time { 20.min * params.readnumber.values().size() * task.attempt }
    maxRetries 3
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }
    publishDir(
        path: "${params.output_dir_br_input}", 
        mode: 'symlink'
    )
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
    awk '{OFS="\\n"; print \$1 " " \$2,\$3,\$4,\$5}' \
    > ${params.exp_name}_combined_barcoded.fq
    """
}