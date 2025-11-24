#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */
// General
params.exp_name = "benchmark_hg002_chr21"

// Variant preprocessing
// params.chromosomes = ["21", "22"]
params.chromosome = "21"
params.link_ref_base = "https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome."//${params.chr}.fa.gz"
params.link_vcf = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
params.link_vcf_idx = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"

// Preprocessing
// software_dir needs to be defined in nextflow.config
// Nanosim
params.abs_path_to_nanosim_repo = "${params.software_dir}/NanoSim"
params.model = "${params.abs_path_to_nanosim_repo}/pre-trained_models/human_giab_hg002_sub1M_kitv14_dorado_v3.2.1/training"
params.readnumber = 1000000

// TOML writing
params.mu = 400
params.dump_time = 35000000
params.accept_unmapped = "true"

// Analysis
params.pores = 512
params.analyse_script_path = "${projectDir}/modules/analyse/scripts"
params.conda = "${projectDir}/modules/analyse/envs"


include {VARIANT_INPUT} from './modules/variant_input/boss_variant_input.nf'
include {PREPROCESS_BOSS} from './modules/preprocess/boss_preprocess.nf'
include {SEQUENCE_BOSS;SEQUENCE_PROFILE_BOSS} from './modules/sequence/boss_sequence.nf'
include {ANALYSE_BOSS;BENCHMARK_VCF} from './modules/analyse/boss_analyse.nf'

workflow  {
    // Get genome with variants from benchmark vcf set
    input_ref = Channel.of(params.link_ref_base)
    input_vcf = Channel.of(params.link_vcf)
    input_vcf_idx = Channel.of(params.link_vcf_idx)

    var_output = VARIANT_INPUT(input_ref, input_vcf, input_vcf_idx)

    ref_genome = var_output.ref
    // benchmark_ground_truth = var_output.subset_vcf -- will use later, once I know how to integrate epi2me pipeline
    ind_genome = var_output.ind_genome

    // Preprocess genome for sequence simulation
    preprocessed = PREPROCESS_BOSS(ind_genome, ref_genome)

    // Run sequence simulation
    sequenced = SEQUENCE_PROFILE_BOSS(preprocessed)

    // Analyse seq output
    ANALYSE_BOSS(sequenced.reads_dir, sequenced.l, ref_genome)
}