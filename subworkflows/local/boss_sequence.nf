#!/usr/bin/env nextflow

include {runBRSim} from "../../modules/local/runBRSim.nf"
include {TimeProfilerunBRSim} from "../../modules/local/TimeProfilerunBRSim.nf"
include {MemProfilerunBRSim} from "../../modules/local/MemProfilerunBRSim.nf"

workflow SEQUENCE_BOSS{
    take:
        input_toml
    main:
        seq = runBRSim(input_toml)
    emit:
        dir = seq.dir
        reads_dir = seq.reads_dir
        l = seq.log
}

workflow SEQUENCE_PROFILE_BOSS{
    take:
        input_toml 
    main:
        seq = TimeProfilerunBRSim(input_toml)
        MemProfilerunBRSim(input_toml)
    emit:
        dir = seq.dir
        reads_dir = seq.reads_dir
        l = seq.log
        lprofile = seq.lprofile

}
