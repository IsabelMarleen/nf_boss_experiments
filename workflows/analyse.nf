include {ANALYSE_BOSS} from '../subworkflows/local/boss_analyse.nf'

workflow  ANALYSE{
    // Input reads
    input_reads = channel.fromPath( "${params.reads}/*.fa" )
    // Sequence log
    seq_log = channel.of("${params.log}")
    // Index reference
    ref_input = channel.fromPath("${params.ref}")

    ANALYSE_BOSS(input_reads, seq_log, ref_input)
}