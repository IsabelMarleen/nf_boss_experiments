#!/usr/bin/env nextflow

/*
 * Pipeline parameters should be included in nextflow.config or nf-params.json
 */


include {PREPROCESS_BOSS} from './subworkflows/local/boss_preprocess.nf'
include {SEQUENCE_BOSS;SEQUENCE_PROFILE_BOSS} from './subworkflows/local/boss_sequence.nf'
include { SIM } from './workflows/simulate_fragments.nf'
include { TRIOS } from './workflows/trios.nf'
include { ANALYSE } from './workflows/analyse.nf'

workflow  {
    main:
        // Determine which workflow to run based on param

        if (params.type == "trios") {  
            TRIOS()
        } else if (params.type == "analyse") {
            ANALYSE()
        } else if (params.type == "simulate") {
            SIM()
        }
}
