#!/usr/bin/env nextflow

// Run using 'nextflow run boss_sequence.nf -entry SEQUENCE_BOSS -with-conda -resume -c ../nextflow.config'
/*
 * Pipeline parameters in nextflow.config
 */


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

process TimeProfilerunBRSim {
    clusterOptions '--mem=64G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=12:00:00'
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

process MemProfilerunBRSim {
    clusterOptions '--mem=64G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=12:00:00'
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

workflow SEQUENCE_BOSS{
    take:
        input_toml
    main:
        seq = runBRSim(input_toml)
    emit:
        dir = seq.out[0]
        reads_dir = seq.out[1]
        log = seq.out[2]
}

workflow SEQUENCE_PROFILE_BOSS{
    take:
        input_toml
    main:
        seq = TimeProfilerunBRSim(input_toml)
        MemProfilerunBRSim(input_toml)
    emit:
        dir = seq.dir
        reads_dir = seq.reads_dir
        l = seq.log
        cprofile = seq.cprofile
        cprofile_stdout = seq.cprofile_stdout

}

workflow  {
    input_toml = channel.of(params.toml)
    SEQUENCE_PROFILE_BOSS(input_toml)
}