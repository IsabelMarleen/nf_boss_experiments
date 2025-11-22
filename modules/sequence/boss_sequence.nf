#!/usr/bin/env nextflow

// Run using 'nextflow run boss_sequence.nf -entry SEQUENCE_BOSS -with-conda -resume -c ../nextflow.config'
/*
 * Pipeline parameters
 */
params.exp_name = "benchmark_hg002_chr21"
params.chr = "21"

params.fasta = "/nfs/research/goldman/ipoetzsch/trio_exp2/GCA_009914755-chromosome-21.fasta"
params.ref_link = "/nfs/research/goldman/ipoetzsch/trio_exp2/GCA_009914755-chromosome-21.fasta"

params.abs_path_to_nanosim_repo = "/hps/software/users/goldman/ipoetzsch/NanoSim"
params.model = "${params.abs_path_to_nanosim_repo}/pre-trained_models/human_giab_hg002_sub1M_kitv14_dorado_v3.2.1/training"
params.readnumber = 1000000

params.abspath_to_boss_runs_repo = "/hps/software/users/goldman/ipoetzsch/BOSS-RUNS2"
params.mu = 400
params.toml = "/hps/nobackup/goldman/ipoetzsch/benchmark_hg002_chr21/br_input/static_benchmark_hg002_chr21.toml"

params.base_out_dir = "/hps/software/users/goldman/ipoetzsch/${params.exp_name}"
params.output_dir_data = "${params.base_out_dir}/data"
params.output_dir_nanosim = "${params.base_out_dir}/nanosim_out"
params.output_dir_br_input = "${params.base_out_dir}/br_input"
params.output_dir_br_output = "${params.base_out_dir}/br_output"

// Run br sim
process runBRSim {
    clusterOptions '--mem=32G --nodes=4 --cpus-per-task=32 --ntasks=1 --time=03:00:00'
    publishDir "${params.output_dir_br_output}", mode: 'symlink'
    conda '/hps/software/users/goldman/ipoetzsch/conda/envs/boss_profile'
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
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=03:00:00'
    publishDir "${params.output_dir_br_output}", mode: 'symlink'
    conda '/hps/software/users/goldman/ipoetzsch/conda/envs/boss_profile'
    cache false
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
    # python -m kernprof -lv -p -m boss --toml ${toml} > l_profiler.stdout
    # source /hps/software/users/goldman/ipoetzsch/conda/etc/profile.d/conda.sh && source /hps/software/users/goldman/ipoetzsch/conda/etc/profile.d/mamba.sh
    # conda activate boss_profile
    # kernprof -lv -m boss.BOSS --toml /hps/nobackup/goldman/ipoetzsch/benchmark_hg002_chr21/br_input/benchmark_hg002_chr21.toml
    """
}


workflow SEQUENCE_BOSS{
    input_toml = channel.of(params.toml)
    // runBRSim(input_toml)
    TimeProfilerunBRSim(input_toml)
}