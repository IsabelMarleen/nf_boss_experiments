#!/usr/bin/env nextflow

include {getRef} from "../../modules/local/getRef.nf"
include {combineRef} from "../../modules/local/combineRef.nf"
include {getVCF} from "../../modules/local/getVCF.nf"
include {subsetVCF} from "../../modules/local/subsetVCF.nf"
include {getConsensus} from "../../modules/local/getConsensus.nf"

workflow VARIANT_INPUT{
    take:
        input_ref
        input_vcf_idx
    main:
        // Download reference file(s) for each chromosome and benchmark vcfs + index 
        chrs = Channel.from(params.chromosomes)
        downloaded_ref = getRef(input_ref.first(), chrs)

        // Combine reference files into one file if necessary
        if (params.chromosomes.size() > 1){
            ref = combineRef(downloaded_ref.collect())
        }
        else{
            ref = downloaded_ref
        }

        // Download benchmark vcf set and index
        downloaded_vcf_idx = getVCF(input_vcf_idx)

        // Subset VCF to just the chromosomes that we are sequencing and index subsetted vcf
        subset_vcf = subsetVCF(downloaded_vcf_idx)
        // Create a consensus sequence for each 'individual'
        ind_genome = getConsensus(subset_vcf, ref.first())
        
    emit:
        ref
        subsetVCF
        ind_genome
}
