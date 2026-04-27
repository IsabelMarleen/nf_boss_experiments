process TimeProfilerunBRSim {
    input: 
        path toml
    output: 
        path("out_*"), emit: dir
        path("00_reads/*.fa"), emit: reads_dir
        path("logs/*_boss.log"), emit: log
        path "*_cprofile", emit: cprofile
        path "*_cprofile.stdout", emit: cprofile_stdout

    script:
    """
    export PROFILE_RUN=True
    python -m cProfile -o "${params.exp_name}_cprofile" -m boss.BOSS --toml ${toml} > "${params.exp_name}_cprofile.stdout"
    """
}