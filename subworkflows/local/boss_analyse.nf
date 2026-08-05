#!/usr/bin/env nextflow

include {visualise_simulation} from "../../modules/local/visualise_simulation.nf"

workflow ANALYSE_BOSS{
    take:
        cov
        unblocks
        analysed_log
        otu
        benchmark
    main:
        // Create plots
        plots = visualise_simulation(cov, unblocks, analysed_log, otu, benchmark)
    emit:
        coverage = cov
        unblocks = unblocks
        plots = plots
}

