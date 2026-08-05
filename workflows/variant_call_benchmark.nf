#!/usr/bin/env nextflow

include { BENCHMARK_HUMAN_VCF } from '../subworkflows/local/boss_human_variantcall.nf'

workflow  VARIANT_CALL_BENCHMARK{
    // input merged bam/bai
    bam_bai = channel.fromPath("${params.mapped_dir}/*")
            .map { tuple( it.simpleName, it ) }
            .groupTuple(sort:{a,b -> if (a.extension == 'bam') {
                                return -1;
                            } else if (b.extension == 'bam') {
                                return 1;
                            } else {
                                return 0;
                            }})
            .map{_it, bam_bai -> tuple([id: bam_bai[0].simpleName], bam_bai[0], bam_bai[1])}
            .filter{_meta, bam, bai -> bam.isFile() & bai.isFile()}
        
    // Reference
    ref_input = channel.fromPath("${params.ref}")

    // Truth vcf
    truth_vcf = channel.of(params.truth_vcf).flatten()
    BENCHMARK_HUMAN_VCF(bam_bai, ref_input, truth_vcf)
}