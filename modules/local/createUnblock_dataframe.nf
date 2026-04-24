process createUnblock_dataframe {
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