#!/usr/bin/env nextflow

// Run using 'nextflow run boss_sequence.nf -entry SEQUENCE_BOSS -with-conda -resume -c ../nextflow.config'
/*
 * Pipeline parameters
 */
params.exp_name = "benchmark_hg002_chr21"

params.abspath_to_boss_runs_repo = "${params.software_dir}/BOSS-RUNS2"

params.base_in_dir = "${params.base_out_dir}/${params.exp_name}"
params.seq_br_output = "${params.base_in_dir}/sequence/br_output"

params.br_input = "${params.base_in_dir}/preprocess/br_input"
// params.toml = "${params.br_input}/static_benchmark_hg002_chr21.toml"
params.toml = "/hps/nobackup/goldman/ipoetzsch/boss_experiments/work/a0/c4632f88faa413f3b321fe0c40c1ce/benchmark_hg002_chr21.toml"


// Run br sim
process runBRSim {
    clusterOptions '--mem=64G --nodes=4 --cpus-per-task=32 --ntasks=1 --time=12:00:00'
    publishDir "${params.seq_br_output}", mode: 'symlink'
    conda "${params.conda_base_dir}/boss_profile"
    input: 
        path toml
    output: 
        path("out_*"), emit: dir
        path("00_reads"), emit: reads_dir
        path("*_boss.log"), emit: log

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
        path("00_reads"), emit: reads_dir
        path("*_boss.log"), emit: log
        path "*_cprofile", emit: cprofile
        path "*_cprofile.stdout", emit: cprofile_stdout
        // path("profile_output.txt"), emit:profile_output
        // path "profile_output_*.txt", emit:profile_time
        // path "profile_output.lprof", emit:profile_lprof
        // path "l_profiler.stdout", emit:stdout

    script:
    """
    PROFILE_RUN=True
    python -m cProfile -o "${params.exp_name}_cprofile" -m boss.BOSS --toml ${toml} > "${params.exp_name}_cprofile.stdout"
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
    emit:
        dir = seq.dir
        reads_dir = seq.reads_dir
        l = seq.log
        cprofile = seq.cprofile
        cprofile_stdout = seq.cprofile_stdout

}

workflow  {
    input_toml = channel.of(params.toml)
    // TimeProfilerunBRSim(input_toml)
    SEQUENCE_PROFILE_BOSS(input_toml)
}