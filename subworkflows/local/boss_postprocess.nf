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
include { BIOBAMBAM_BAMMERGE } from '../../modules/nf-core/biobambam/bammerge/main'

workflow POSTPROCESS_BOSS{
    take:
        input_reads
        seq_log
        ref_input
    main:
        // Input reads            
        input_reads
            .flatten()
            .filter { it.simpleName.split('_')[1] != "0" }
            .map { tuple( it.simpleName, it ) }
            .groupTuple()
            .set { input_reads_collection }
        // Index reference
        processed_ref = indexPaf(ref_input)
        ref_idx = processed_ref.ref_idx
        ref_unzipped = processed_ref.ref_unzipped

        // Map paf
        first_map = mapPaf(input_reads.flatten(), ref_idx.collect())
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
        sep_target = separateTarget(first_map_tuple, otu.collect())

        // Map sam
        mapped_sam = mapSam(sep_target.flatten(), ref_unzipped.collect(), otu.collect())

        bam_to_merge = mapped_sam
            .map { pair ->
                def bam  = pair[0]
                def base = bam.baseName

                def parts = base.split('_') as List

                def type = null
                def idx = null
                def type2 = null
                def type3 = null
                if (parts.size() == 3) {
                // type_index_type2
                type  = parts[0]
                idx   = parts[1] as int
                type2 = parts[2]
                if( idx == 0 ) return null
                }
                else if (parts.size() == 4) {
                // type_index_type2_type3
                type  = parts[0]
                idx   = parts[1] as int
                type2 = parts[2]
                type3 = parts[3]
                if( idx == 0 ) return null
                }
                else {
                throw new IllegalArgumentException("Unexpected BAM filename format: ${base}")
                }

                def key = (type3 == null) ? "${type}_${type2}" : "${type}_${type2}_${type3}"

                // Important: make the grouped value a single element we can sort
                [ key, [idx, bam] ]
            }
            .groupTuple()  // emits: [ key, listOfValues ], where each value is [idx, bam]
            .flatMap { _key, values ->
                def sorted = values.sort { it[0] }        // sort by idx
                def bams = sorted.collect { it[1] }      // extract bam paths

                // cumulative prefixes
                (0..<bams.size()).collect { i -> tuple([id:bams[i].simpleName+'_merged'], bams[0..i]) }
        }

        BIOBAMBAM_BAMMERGE(bam_to_merge)
        
        merged_bam = BIOBAMBAM_BAMMERGE.out.bam
            .combine(BIOBAMBAM_BAMMERGE.out.bam_index, by:0)

        // Pileup
        pile = pileup(merged_bam, otu.collect())

        // Record unblocks
        unblocks_ind = recordUnblocks(sep_target.flatten())

        // Create unblocks dataframe
        unblocks = createUnblock_dataframe(unblocks_ind.collect())

        // Create coverage dataframe
        cov = create_coverage_dataframe(pile.csv.collect())

        // Analyse log
        analysed_log = analyse_log(seq_log)

    emit:
        merged_bam = merged_bam
        coverage = cov
        unblocks = unblocks
        l = analysed_log
        otu = otu
        ref_unzipped
}

