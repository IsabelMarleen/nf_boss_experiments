process create_coverage_dataframe { 
    input:
        path csv
    output:
        path 'coverage.csv', emit: coverage

    script: 
    if (params.barcodes == null){
        bc_string = ""
    }
    else{
        bc_string = "bc,"
    }
    """
    mkdir -p results/$params.exp_name &&
    cat <(echo cond,time,otu,${bc_string}mean_coverage,low_coverage_prop,evenness) $csv > coverage.csv
    """
} 