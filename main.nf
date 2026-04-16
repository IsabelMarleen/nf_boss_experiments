#!/usr/bin/env nextflow

/*
 * Pipeline parameters should be included in nextflow.config or nf-params.json
 */


include {VARIANT_INPUT} from './workflows/boss_variant_input.nf'
include {SIMULATE_FRAGMENTS_BOSS} from './workflows/boss_simulate_fragments.nf'
include {PREPROCESS_BOSS} from './workflows/boss_preprocess.nf'
include {SEQUENCE_BOSS;SEQUENCE_PROFILE_BOSS} from './workflows/boss_sequence.nf'
include {ANALYSE_BOSS} from './workflows/boss_analyse.nf'

workflow  {
    // Get genome with variants from benchmark vcf set
    input_ref = Channel.of(params.link_ref)
    input_vcf = Channel.of(params.link_vcf).flatMap()

    var_output = VARIANT_INPUT(input_ref, input_vcf)

    ref_genome = var_output.ref
    // benchmark_ground_truth = var_output.subset_vcf -- will use later, once I know how to integrate epi2me pipeline
    ind_genome = var_output.ind_genome

    // Preprocess genome for sequence simulation
    sim_fragments = SIMULATE_FRAGMENTS_BOSS(ind_genome)

    // Preprocess genome for sequence simulation
    preprocessed = PREPROCESS_BOSS(sim_fragments, ref_genome)

    // Run sequence simulation
    sequenced = SEQUENCE_PROFILE_BOSS(preprocessed)

    // Analyse seq output
    ANALYSE_BOSS(sequenced.reads_dir, sequenced.l, ref_genome)
}