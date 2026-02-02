#!/usr/bin/env nextflow

/*
 * Pipeline parameters should be included in nextflow.config or nf-params.json
 */


include {VARIANT_INPUT} from './modules/boss_variant_input.nf'
include {PREPROCESS_BOSS} from './modules/boss_preprocess.nf'
include {SEQUENCE_BOSS;SEQUENCE_PROFILE_BOSS} from './modules/boss_sequence.nf'
include {ANALYSE_BOSS;BENCHMARK_VCF} from './modules/boss_analyse.nf'

workflow  {
    // Get genome with variants from benchmark vcf set
    input_ref = Channel.of(params.link_ref)
    input_vcf = Channel.of(params.link_vcf).flatMap()

    var_output = VARIANT_INPUT(input_ref, input_vcf)

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