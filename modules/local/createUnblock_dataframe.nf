process createUnblock_dataframe {
    clusterOptions '--mem=12G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_debug}/results", mode: 'symlink'
    input:
        path csv
    output:
        path 'unblocks.csv', emit:unblocks
    script:
    if (params.barcodes == null){
        bc_string = ""
    }
    else{
        bc_string = "bc,"
    }
    """
    cat <(echo cond,time,otu,${bc_string}total,base_total,unb,unb_ratio) $csv > unblocks.csv
    """
}