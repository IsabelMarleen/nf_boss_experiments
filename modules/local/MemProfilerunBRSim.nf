process MemProfilerunBRSim {
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
        path "*.bin", emit: memprofile

    script:
    """
    export PROFILE_RUN=True
    memray run --follow-fork --aggregate -m boss.BOSS --toml ${toml}
    """
}