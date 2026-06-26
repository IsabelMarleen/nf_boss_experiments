#!/usr/bin/env nextflow

include {indexPaf} from "../../modules/local/indexPaf.nf"
include {mapPaf} from "../../modules/local/mapPaf.nf"
include {processOtu} from "../../modules/local/processOtu.nf"
include {separateTarget} from "../../modules/local/separateTarget.nf"
include {mapSam} from "../../modules/local/mapSam.nf"
include {pileup} from "../../modules/local/pileup.nf"
include {recordUnblocks} from "../../modules/local/recordUnblocks.nf"
include {createUnblock_dataframe} from "../../modules/local/createUnblock_dataframe.nf"
include {create_coverage_dataframe} from "../../modules/local/create_coverage_dataframe.nf"
include {analyse_log} from "../../modules/local/analyse_log.nf"
include {visualise_simulation} from "../../modules/local/visualise_simulation.nf"

workflow ANALYSE_BOSS{
    take:
        input_reads
        seq_log
        ref_input
    main:
        // Input reads            
        input_reads
            .flatten()
            .map { tuple( it.simpleName, it ) }
            .groupTuple()
            .set { input_reads_collection }
        // Index reference
        processed_ref = indexPaf(ref_input)
        ref_idx = processed_ref.ref_idx
        ref_unzipped = processed_ref.ref_unzipped

        // Map paf
        first_map = mapPaf(input_reads.flatten(), ref_idx.first())
        first_map
                .map { tuple( it.simpleName, it ) }
                .combine( input_reads_collection, by: 0 )
                .transpose( by: 2 )
                .map { _case_id, paf, fa -> tuple( paf, fa ) }
                .set{first_map_tuple}

        // Create otu_info.py file
        otu_names = channel.of(params.otu_clean_names)
        otu_sizes = channel.of(params.otu_genome_sizes)
        otu = processOtu(otu_names, otu_sizes)
        // Separate target
        sep_target = separateTarget(first_map_tuple, otu.first())

        // Map sam
        mapped_sam = mapSam(sep_target.flatten(), ref_unzipped.first(), otu.first())

        // If barcoded experiment, add barcode info to each fq and combine
        if (params.barcodes != null){
            pile_bam = mapped_sam
                .filter { it[0].simpleName.split('_')[1] != "0" }
                .toSortedList{it -> it[0].simpleName.split('_')[1].toInteger()}
                .flatMap()
                .map{tuple(it[0].simpleName.split('_')[0,2,3], it[0])}
                .groupTuple()

            pile_bai = mapped_sam
                .filter { it[0].simpleName.split('_')[1] != "0" }
                .toSortedList{it -> it[0].simpleName.split('_')[1].toInteger()}
                .flatMap()
                .map{tuple(it[0].simpleName.split('_')[0,2,3], it[1])}
                .groupTuple()
        }
        else{
            pile_bam = mapped_sam
                .filter { it[0].simpleName.split('_')[1] != "0" }
                .toSortedList{it -> it[0].simpleName.split('_')[1].toInteger()}
                .flatMap()
                .map{tuple(it[0].simpleName.split('_')[0,2], it[0])}
                .groupTuple()   
            pile_bai = mapped_sam
                .filter { it[0].simpleName.split('_')[1] != "0" }
                .toSortedList{it -> it[0].simpleName.split('_')[1].toInteger()}
                .flatMap()
                .map{tuple(it[0].simpleName.split('_')[0,2], it[1])}
                .groupTuple()  
        }
        // Pileup
        pile = pileup(pile_bam, pile_bai, otu.first())

        // Record unblocks
        unblocks_ind = recordUnblocks(sep_target.flatten())

        // Create unblocks dataframe
        unblocks = createUnblock_dataframe(unblocks_ind.collect())

        // Create coverage dataframe
        cov = create_coverage_dataframe(pile.csv.collect())

        // Analyse log
        analysed_log = analyse_log(seq_log)

        // Create plots
        plots = visualise_simulation(cov, unblocks, analysed_log, otu)
    emit:
        coverage = cov
        unblocks = unblocks
        plots = plots
}

