#!/usr/bin/env nextflow

include {VARIANT_INPUT} from './modules/variant_input/boss_variant_input.nf'
include {PREPROCESS_BOSS} from './modules/preprocess/boss_preprocess.nf'
include {SEQUENCE_BOSS} from './modules/sequence/boss_sequence.nf'
include {ANALYSE_BOSS} from './modules/analyse/boss_analyse.nf'