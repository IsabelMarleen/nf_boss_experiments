process MemProfilerunBRSim {
    input: 
        path toml
    output: 
        path("out_*"), emit: dir
        path("00_reads/*.fa"), emit: reads_dir
        path("logs/*_boss.log"), emit: log
        path "*.bin", emit: memprofile

    script:
    """
    export PROFILE_RUN=True
    memray run --follow-fork --aggregate -m boss.BOSS --toml ${toml}
    """
}