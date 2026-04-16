#!/usr/bin/env nextflow

// Run using 'nextflow run boss_preprocess.nf -entry PREPROCESS_BOSS -with-conda -resume -c ../nextflow.config'

/*
 * Pipeline parameters in ../nextflow.config
 */


include {truncateFq} from "../modules/local/truncateFq.nf"
include {scan_offsetsFq} from "../modules/local/scan_offsetsFq.nf"
include {mapPaf} from "../modules/local/mapPaf.nf"
include {mapPaf_trunc} from "../modules/local/mapPaf_trunc.nf"
include {scan_offsets_Paf} from "../modules/local/scan_offsets_Paf.nf"
include {writeToml} from "../modules/local/writeToml.nf"

workflow PREPROCESS_BOSS{
    take:
        aligned_fastq
        ref
    main:
        // Prepare input for BRsim
        // Truncate fq Step 1 of prepare input for BRsim
        trunc_fq = truncateFq(aligned_fastq)

        // Scan fq offsets Step 2 of prepare input for BRsim
        boss_sampler = Channel.fromPath(params.abs_path_to_boss_repo+"/sampler.py")
        fq_combined = trunc_fq.full_fq.mix(trunc_fq.trunc_fq)
        fq_offsets = scan_offsetsFq(fq_combined, boss_sampler)
        // TODO: Combine input genomes into one file using cat

        // Map paf Step 3 of prepare input for BRsim
        mpaf = mapPaf(aligned_fastq, ref)

        // Map paf trunc Step 4 of prepare input for BRsim
        mpaf_trunc = mapPaf_trunc(trunc_fq.trunc_fq, ref)

        // Scan paf offsets Step 5 of prepare input for BRsim
        paf_offsets = scan_offsets_Paf(mpaf.mappings, mpaf_trunc.mappings_trunc)

        // Create toml file based on input files above
        toml = writeToml(trunc_fq.full_fq, mpaf, mpaf_trunc, trunc_fq.trunc_fq, fq_offsets.fq_offsets.collect(), 
                        paf_offsets.mappings_offsets, paf_offsets.mappings_offsets_trunc, ref)

        // runBRSim(toml)
    emit:
        toml
}


workflow PREPROCESS_From_seq{
    take:
        full_fq
        ref
    main:
        // Prepare input for BRsim
        // Truncate fq Step 1 of prepare input for BRsim
        trunc_fq = truncateFq(full_fq)

        // Scan fq offsets Step 2 of prepare input for BRsim
        boss_sampler = Channel.fromPath(params.abs_path_to_boss_repo+"/sampler.py")
        fq_combined = full_fq.mix(trunc_fq.trunc_fq)
        fq_offsets = scan_offsetsFq(fq_combined, boss_sampler)

        // Map paf Step 3 of prepare input for BRsim
        mpaf = mapPaf(full_fq, ref)

        // Map paf trunc Step 4 of prepare input for BRsim
        mpaf_trunc = mapPaf_trunc(trunc_fq.trunc_fq, ref)

        // Scan paf offsets Step 5 of prepare input for BRsim
        paf_offsets = scan_offsets_Paf(mpaf.mappings, mpaf_trunc.mappings_trunc)

        // Create toml file based on input files above
        toml = writeToml(full_fq, mpaf, mpaf_trunc, trunc_fq.trunc_fq, fq_offsets.fq_offsets.collect(), 
                        paf_offsets.mappings_offsets, paf_offsets.mappings_offsets_trunc, ref)

        // runBRSim(toml)
    emit:
        toml
}

// workflow  PREPROCESS_WITH_SEQ{
//     input_genome = Channel.of(params.genome)
//     reference = Channel.of(params.ref)
//     // Run preprocessing
//     toml = PREPROCESS_BOSS(input_genome, reference) 

// }

// workflow  {
//     // Download reference file(s) for each chromosome and benchmark vcfs + index 
//     input_genome = Channel.of(params.genome)
//     reference = Channel.of(params.ref)

//     PREPROCESS_BOSS(input_genome, reference)
// }

workflow  {
    // Download reference file(s) for each chromosome and benchmark vcfs + index 
    reference = Channel.fromPath("${params.abs_path_to_boss_repo}/data/zymo.fa")
    full_fq = Channel.fromPath("${params.software_dir}/boss_exp/barcoded_input_data/FAT91932_pass_e7bf7751_f43c451e_4_bc.fastq")

    PREPROCESS_From_seq(full_fq, reference)
}