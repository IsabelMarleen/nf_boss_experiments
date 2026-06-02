#!/usr/bin/env nextflow

include {truncateFq} from "../../modules/local/truncateFq.nf"
include {scan_offsetsFq} from "../../modules/local/scan_offsetsFq.nf"
include {mapPaf} from "../../modules/local/mapPaf.nf"
include {mapPaf_trunc} from "../../modules/local/mapPaf_trunc.nf"
include {scan_offsets_Paf} from "../../modules/local/scan_offsets_Paf.nf"
include {writeToml} from "../../modules/local/writeToml.nf"

workflow PREPROCESS_BOSS{
    take:
        aligned_fastq
        ref
    main:
        // Prepare input for BRsim
        // Truncate fq Step 1 of prepare input for BRsim
        trunc_fq = truncateFq(aligned_fastq)

        // Scan fq offsets Step 2 of prepare input for BRsim
        fq_combined = trunc_fq.full_fq.mix(trunc_fq.trunc_fq)
        fq_offsets = scan_offsetsFq(fq_combined)

        // Map paf Step 3 of prepare input for BRsim
        mpaf = mapPaf(aligned_fastq, ref)

        // Map paf trunc Step 4 of prepare input for BRsim
        mpaf_trunc = mapPaf_trunc(trunc_fq.trunc_fq, ref)

        // Scan paf offsets Step 5 of prepare input for BRsim
        paf_offsets = scan_offsets_Paf(mpaf.mappings, mpaf_trunc.mappings_trunc)

        // Create toml file based on input files above
        toml = writeToml(trunc_fq.full_fq, mpaf, mpaf_trunc, trunc_fq.trunc_fq, fq_offsets.fq_offsets.collect(), 
                        paf_offsets.mappings_offsets, paf_offsets.mappings_offsets_trunc, ref)

    emit:
        toml
}
