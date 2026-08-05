#!/usr/bin/env nextflow


include {simNanofastq} from "../../modules/local/simNanofastq.nf"
include {addBarcodes} from "../../modules/local/addBarcodes.nf"
include {combineBarcodedFq} from "../../modules/local/combineBarcodedFq.nf"
include {combineLargeBarcodedFq} from "../../modules/local/combine_large_bc_fq.nf"

workflow SIMULATE_FRAGMENTS_BOSS{
    take:
        input_genome
    main:
        // Simulate reads using Nanosim
        nanosim = simNanofastq(input_genome)

        // If barcoded experiment, add barcode info to each fq and combine
        if (params.barcodes != null){
            barcoding_needed = channel.of(true)
            barcoded = addBarcodes(nanosim.aligned_fastq, barcoding_needed.collect())
            if ( params.readnumber.values().sum() < 9000000){
                aligned_fastq = combineBarcodedFq(barcoded.collect())
            }
            else{
                aligned_fastq = combineLargeBarcodedFq(barcoded.collect())
            }
                
        }
        else{
            aligned_fastq = nanosim.aligned_fastq
        }

       
    emit:
        aligned_fastq
}
