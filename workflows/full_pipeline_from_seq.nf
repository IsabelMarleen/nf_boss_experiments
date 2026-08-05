#!/usr/bin/env nextflow

include {PREPROCESS_BOSS} from '../subworkflows/local/boss_preprocess.nf'
include {SEQUENCE_PROFILE_BOSS; SEQUENCE_BOSS} from '../subworkflows/local/boss_sequence.nf'
include {ANALYSE_BOSS} from '../subworkflows/local/boss_analyse.nf'
include {POSTPROCESS_BOSS} from '../subworkflows/local/boss_postprocess.nf'
include { BENCHMARK_HUMAN_VCF } from '../subworkflows/local/boss_human_variantcall.nf'
include {getMapRef} from '../modules/local/getMapRef.nf'
include { getUnzipComb } from '../modules/local/getUnzipComb.nf'
include { getGunzipReads } from '../modules/local/getGunzipReads.nf'
include { combineLargeBarcodedFq } from '../modules/local/combine_large_bc_fq.nf'
include { combineBarcodedFq } from '../modules/local/combineBarcodedFq.nf'
include { addBarcodes } from '../modules/local/addBarcodes.nf'


workflow FROMSEQ{
    main:
        // Get genome with variants from benchmark vcf set
        input_ref = Channel.of(params.link_ref)
        input_reads = Channel.of(params.reads)

        // Download ref, Might have to combine ref
        ref_genome = getUnzipComb(input_ref.flatten())
        // Download reads, Might have to gunzip reads
        unz_reads = getGunzipReads(input_reads.flatten())

        if (params.reads.size() > 1 ){
            def barcoding_needed_ch = unz_reads.map { f ->
                def need = false
                f.withReader { r ->
                    def firstLine = r.readLine()
                    need = (firstLine == null) || !firstLine.contains('barcode=')
                }
                need
            }
            barcoded = addBarcodes(unz_reads, barcoding_needed_ch)
            reads = combineBarcodedFq(barcoded.collect())
        }
        else {
            reads = unz_reads
        }

        // Preprocess genome for sequence simulation
        preprocessed = PREPROCESS_BOSS(reads, ref_genome)

        // Run sequence simulation
        if (params.profile == true){
            sequenced = SEQUENCE_PROFILE_BOSS(preprocessed)
        }
        else {
            sequenced = SEQUENCE_BOSS(preprocessed)
        }

    processed = POSTPROCESS_BOSS(sequenced.reads_dir, sequenced.l, ref_genome)
    if (params.benchmark){
        print params.benchmark
        truth_vcf = channel.of(params.truth_vcf).flatten()
        benchmark_summary = BENCHMARK_HUMAN_VCF(processed.merged_bam, ref_genome, truth_vcf)
    } else{
        benchmark_summary = []
    }
    ANALYSE_BOSS(processed.coverage, processed.unblocks, processed.l, processed.otu, benchmark_summary)
}