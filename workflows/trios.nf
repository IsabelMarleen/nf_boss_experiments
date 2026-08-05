#!/usr/bin/env nextflow

include {VARIANT_INPUT_BOSS} from '../subworkflows/local/boss_variant_input.nf'
include {SIMULATE_FRAGMENTS_BOSS} from '../subworkflows/local/boss_simulate_fragments.nf'
include {PREPROCESS_BOSS} from '../subworkflows/local/boss_preprocess.nf'
include {SEQUENCE_PROFILE_BOSS} from '../subworkflows/local/boss_sequence.nf'
include {ANALYSE_BOSS} from '../subworkflows/local/boss_analyse.nf'
include {POSTPROCESS_BOSS} from '../subworkflows/local/boss_postprocess.nf'
include { BENCHMARK_HUMAN_VCF } from '../subworkflows/local/boss_human_variantcall.nf'
include {getMapRef} from '../modules/local/getMapRef.nf'

workflow TRIOS{
    main:
        // Get genome with variants from benchmark vcf set
        input_ref = Channel.of(params.link_ref)
        input_vcf = Channel.of(params.link_vcf).flatMap()

        var_output = VARIANT_INPUT_BOSS(input_ref, input_vcf)

        // Check whether a different reference should be used for mapping than the one for variant creation
        if (params.link_ref_map != null){
            map_ref = Channel.of(params.link_ref_map)
            ref_genome = getMapRef(map_ref)
        }
        else{
            ref_genome = var_output.ref
        }
        ind_genome = var_output.ind_genome

        // Preprocess genome for sequence simulation
        sim_fragments = SIMULATE_FRAGMENTS_BOSS(ind_genome)

        // Preprocess genome for sequence simulation
        preprocessed = PREPROCESS_BOSS(sim_fragments, ref_genome)

        // Run sequence simulation
        sequenced = SEQUENCE_PROFILE_BOSS(preprocessed)
        
        // Postprocess sequencing results
        processed = POSTPROCESS_BOSS(sequenced.reads_dir, sequenced.log, ref_genome)

        // If applicable, run variant call benchmark
        if (params.benchmark){
            benchmark_summary = BENCHMARK_HUMAN_VCF(processed.merged_bam, ref_genome, var_output.subsetVCF)
        } else{
            benchmark_summary = []
        }

        // Visualise analysis results
        ANALYSE_BOSS(processed.coverage, processed.unblocks, processed.l, processed.otu, benchmark_summary)
}