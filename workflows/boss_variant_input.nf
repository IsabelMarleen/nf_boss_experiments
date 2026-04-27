#!/usr/bin/env nextflow

// Run using 'nextflow run boss_variant_input.nf -entry VARIANT_INPUT -with-conda -resume -c ../nextflow.config'

/*
 * Pipeline parameters in nextflow.config
 */

include {getRef} from "../modules/local/getRef.nf"
include {combineRef} from "../modules/local/combineRef.nf"
include {getVCF} from "../modules/local/getVCF.nf"
include {subsetVCF} from "../modules/local/subsetVCF.nf"
include {getConsensus} from "../modules/local/getConsensus.nf"

workflow VARIANT_INPUT{
    take:
        input_ref
        input_vcf_idx
    main:
        // Download reference file(s) for each chromosome and benchmark vcfs + index 
        // TODO: To enable getting the consensus genomes for several individuals,
        // I could take the reference processing steps, put them into a separate workflow 
        // and then I can call the remaining part of the workflow for each individual
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
        // subset_vcf.view()
        // Create a consensus sequence for each 'individual'
        ind_genome = getConsensus(subset_vcf, ref.first())
        
    emit:
        ref
        subsetVCF
        ind_genome
}

// workflow VARIANT_INPUT_per_bc{
//     take:
//         ref
//         input_vcf
//     main:
//         // Download benchmark vcf set and index
//         downloaded_vcf_idx = getVCF(input_vcf)


//         // Subset VCF to just the chromosomes that we are sequencing and index subsetted vcf
//         subset_vcf = subsetVCF(downloaded_vcf_idx)

//         // Create a consensus sequence for each 'individual'
//         ind_genome = getConsensus(subset_vcf, ref)
        
//     emit:
//         subsetVCF
//         ind_genome
// } 

// workflow barcoded_VARIANT_INPUT{
//     take:
//         input_ref
//         input_vcf_bc1
//         input_vcf_idx_bc1
//         input_vcf_bc2
//         input_vcf_idx_bc2
//     main:
//         // Download reference file(s) for each chromosome and benchmark vcfs + index 
//         // TODO: To enable getting the consensus genomes for several individuals,
//         // I could take the reference processing steps, put them into a separate workflow 
//         // and then I can call the remaining part of the workflow for each individual
//         chrs = Channel.from(params.chromosomes)
//         downloaded_ref = getRef(input_ref.first(), chrs)

//         // Combine reference files into one file if necessary
//         if (params.chromosomes.size() > 1){
//             ref = combineRef(downloaded_ref.collect())
//         }
//         else{
//             ref = downloaded_ref
//         }

//         // VARIANT_INPUT_per_bc(ref, )
        
//     emit:
//         ref
//         subsetVCF
//         // ind_genome
// } 

workflow {
    input_ref = Channel.of(params.link_ref)
    input_vcf = Channel.of(params.link_vcf).flatMap()

    VARIANT_INPUT(input_ref, input_vcf)
}