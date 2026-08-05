process TimeProfilerunBRSim {
    input: 
        path toml
    output: 
        path("out_*"), emit: dir
        path("00_reads/*.fa"), emit: reads_dir
        path("logs/*_boss.log"), emit: log
        path("*.lprof"), emit: lprofile
    script:
    """
    python -m kernprof -lv -p -m boss.BOSS --toml ${toml}
    """
}