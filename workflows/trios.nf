#!/usr/bin/env nextflow

include {VARIANT_INPUT_BOSS} from '../subworkflows/local/boss_variant_input.nf'
include {SIMULATE_FRAGMENTS_BOSS} from '../subworkflows/local/boss_simulate_fragments.nf'
include {PREPROCESS_BOSS} from '../subworkflows/local/boss_preprocess.nf'
include {SEQUENCE_PROFILE_BOSS} from '../subworkflows/local/boss_sequence.nf'
include {ANALYSE_BOSS} from '../subworkflows/local/boss_analyse.nf'

workflow TRIOS{
    main:
        // Get genome with variants from benchmark vcf set
        input_ref = Channel.of(params.link_ref)
        input_vcf = Channel.of(params.link_vcf).flatMap()

        var_output = VARIANT_INPUT_BOSS(input_ref, input_vcf)

        ref_genome = var_output.ref
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