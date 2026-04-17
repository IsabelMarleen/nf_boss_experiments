// Visualise results
process visualise_simulation {
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:20:00'
    publishDir "${params.output_dir_debug}/results", mode: 'copy'
    conda "${params.conda_envs}/visualisation.yaml"
    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d5f38f55a40143553c791afd3d338ccc13c7de139c1a1d6ca78f3e389c49c89b/data"
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
    visualise_simulation.R \
    --input_cov ${coverage} \
    --input_unb ${unblocks} \
    --dump_time ${params.dump_time} \
    --pores ${params.pores} \
    --analysed_log $analysed_log \
    --output ${params.exp_name}_plots.pdf
    """
}