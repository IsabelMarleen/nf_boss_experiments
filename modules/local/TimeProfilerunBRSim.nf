process TimeProfilerunBRSim {
    clusterOptions '--nodes=1 --cpus-per-task=32 --ntasks=1'
    memory { 64.GB * params.readnumber.values().size() * params.chromosomes.size() * task.attempt }
    time {6.h * params.readnumber.values().size() * params.chromosomes.size() * task.attempt}
    publishDir "${params.seq_br_output}", mode: 'symlink'
    conda "${params.conda_base_dir}/boss_profile"
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