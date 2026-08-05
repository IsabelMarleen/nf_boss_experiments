#!/usr/bin/env nextflow

include {ANALYSE_BOSS} from '../subworkflows/local/boss_analyse.nf'
include {POSTPROCESS_BOSS} from '../subworkflows/local/boss_postprocess.nf'
include { BENCHMARK_HUMAN_VCF } from '../subworkflows/local/boss_human_variantcall.nf'

workflow  ANALYSE{
    // Input reads
    input_reads = channel.fromPath( "${params.reads}/*.fa" )
    // Sequence log
    seq_log = channel.of("${params.log}")
    // Index reference
    ref_input = channel.fromPath("${params.ref}")
    

    processed = POSTPROCESS_BOSS(input_reads, seq_log, ref_input)
    if (params.benchmark){
        truth_vcf = channel.of(params.truth_vcf).flatten()
        benchmark_summary = BENCHMARK_HUMAN_VCF(processed.merged_bam, ref_input, truth_vcf)
    } else{
        benchmark_summary = []
    }
    ANALYSE_BOSS(processed.coverage, processed.unblocks, processed.l, processed.otu, benchmark_summary)

}