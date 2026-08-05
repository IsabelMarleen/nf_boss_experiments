#!/usr/bin/env nextflow

include {getFile} from "../../modules/local/getFile.nf"
include {combineRef} from "../../modules/local/combineRef.nf"
include {getVCF} from "../../modules/local/getVCF.nf"
include {subsetVCF} from "../../modules/local/subsetVCF.nf"
include {getConsensus} from "../../modules/local/getConsensus.nf"

workflow VARIANT_INPUT_BOSS{
    take:
        input_ref
        input_vcf_idx
    main:
        // Download reference file(s) for each chromosome and benchmark vcfs + index 
        chrs = Channel.from(params.chromosomes)

        if (params.link_ref.toString().endsWith("chromosome.")){
            input_ref = chrs
                .map{chr -> "${params.link_ref}${chr}.fa.gz"}
        }

        downloaded_ref = getFile(input_ref)

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
        ind_genome = getConsensus(subset_vcf, ref.collect())
        
    emit:
        ref
        subsetVCF
        ind_genome
}
