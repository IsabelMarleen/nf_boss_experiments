// Run br sim
process runBRSim {
    input: 
        path toml
    output: 
        path("out_*"), emit: dir
        path("00_reads/*.fa"), emit: reads_dir
        path("logs/*_boss.log"), emit: log

    script:
    """
    boss --toml $toml
    """
}