// Run br sim
process runBRSim {
    clusterOptions '--mem=64G --nodes=4 --cpus-per-task=32 --ntasks=1 --time=12:00:00'
    publishDir "${params.seq_br_output}", mode: 'symlink'
    conda "${params.conda_base_dir}/boss_profile"
    input: 
        path toml
    output: 
        path("out_*"), emit: dir
        path("00_reads/*.fa"), emit: reads_dir
        path("logs/*_boss.log"), emit: log

    script:
    """
    boss --toml ${params.toml}
    """
}