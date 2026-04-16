// Visualise results
process visualise_simulation {
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:20:00'
    publishDir "${params.output_dir_debug}/results", mode: 'copy'
    conda "${params.conda_envs}/visualisation.yaml"
    cache false
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
    ${params.script_dir}visualise_simulation.R \
    --input_cov ${coverage} \
    --input_unb ${unblocks} \
    --dump_time ${params.dump_time} \
    --pores ${params.pores} \
    --analysed_log $analysed_log \
    --output ${params.exp_name}_plots.pdf
    """
}