// Visualise results
process visualise_simulation {
    input: 
        path coverage
        path unblocks
        path analysed_log
        path otu
    output: 
        path "*.pdf", emit: result_plot

    script:
    """
    #now=date +'%Y/%m/%d'
    which Rscript
    visualise_simulation.R \
    --input_cov ${coverage} \
    --input_unb ${unblocks} \
    --dump_time ${params.dump_time} \
    --pores ${params.pores} \
    --analysed_log $analysed_log \
    --output ${params.exp_name}_plots.pdf
    """
}