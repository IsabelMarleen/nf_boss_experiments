// Visualise results
process visualise_simulation {
    input: 
        path coverage
        path unblocks
        path analysed_log
        path otu
        path vcf_summary
    output: 
        path "*.pdf", emit: result_plot

    script:
    chr_num = params.otu_clean_names.size()
    """
    #now=date +'%Y/%m/%d'
    which Rscript
    if [[ $chr_num > 6 ]]
    then
        visualise_all_chr.R \
        --input_cov ${coverage} \
        --input_unb ${unblocks} \
        --dump_time ${params.dump_time} \
        --pores ${params.pores} \
        --analysed_log $analysed_log \
        --output ${params.exp_name}_plots.pdf
    else
        visualise_simulation.R \
        --input_cov ${coverage} \
        --input_unb ${unblocks} \
        --dump_time ${params.dump_time} \
        --pores ${params.pores} \
        --analysed_log $analysed_log \
        --output ${params.exp_name}_plots.pdf \
        --output_vcf ${params.exp_name}_vcf_plots.pdf \
        --vcf_summary ${vcf_summary}
    fi
    """
}