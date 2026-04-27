#!/usr/bin/env nextflow


/*
 * Pipeline parameters in ../nextflow.config
 */
include {simNanofastq} from "../modules/local/simNanofastq.nf"
include {addBarcodes} from "../modules/local/addBarcodes.nf"
include {combineBarcodedFq} from "../modules/local/combineBarcodedFq.nf"

workflow SIMULATE_FRAGMENTS_BOSS{
    take:
        input_genome
    main:
        // Simulate reads using Nanosim
        nanosim = simNanofastq(input_genome)

        // If barcoded experiment, add barcode info to each fq and combine
        if (params.barcodes != null){
            barcoded = addBarcodes(nanosim.aligned_fastq)
            aligned_fastq = combineBarcodedFq(barcoded.collect())
        }
        else{
            aligned_fastq = nanosim.aligned_fastq
        }

       
    emit:
        aligned_fastq
}
