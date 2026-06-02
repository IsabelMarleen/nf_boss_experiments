include {VARIANT_INPUT} from '../subworkflows/local/boss_variant_input.nf'
include {SIMULATE_FRAGMENTS_BOSS} from '../subworkflows/local/boss_simulate_fragments.nf'

workflow SIM{
    main:
        // Get genome with variants from benchmark vcf set
        input_ref = Channel.of(params.link_ref)
        input_vcf = Channel.of(params.link_vcf).flatMap()

        var_output = VARIANT_INPUT(input_ref, input_vcf)

        ref_genome = var_output.ref
        ind_genome = var_output.ind_genome

        // Preprocess genome for sequence simulation
        sim_fragments = SIMULATE_FRAGMENTS_BOSS(ind_genome)
}