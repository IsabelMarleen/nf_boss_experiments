#!/usr/bin/env nextflow

// Run using 'nextflow run boss_analyse.nf -with-conda -resume'
// Depends on a nextflow.config file that defines the parameters

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


// process visualise_time {
//     clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:20:00'
//     publishDir "${params.output_dir_debug}/results", mode: 'symlink'
//     conda "${params.conda_base_dir}/boss_profile"
//     input:
//         path time_profile

//     output:
//         path "*.png"

//     script:
//     """
//     gprof2dot -f pstats ${time_profile} | dot -Tpng -o ${time_profile.getSimpleName()}_dot.png
//     """
// }

// process benchmark_variants {
//     tag "${test_vcf.getSimpleName()}"
//     clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:20:00'
//     publishDir "${params.output_dir_debug}/results", mode: 'symlink'
//     input:
//         tuple path(ground_truth_vcf), path(ground_truth_vcf_idx)
//         tuple path(test_vcf), path(test_vcf_idx)
//         path ref

//     output:
//         path "vcf_benchmark_*"
//         path "*.sdf"

//     script:
//     """
//     "${params.rtg_dir_path}/rtg" format -o "${ref.getSimpleName()}.sdf" ${ref}
//     "${params.rtg_dir_path}/rtg" vcfeval -b ${ground_truth_vcf} -c ${test_vcf} -t "${ref.getSimpleName()}.sdf" -o "vcf_benchmark_${test_vcf.getSimpleName()}"
//     """
// }

// process plot_rocs {
//     clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:20:00'
//     publishDir "${params.output_dir_results}", mode: 'symlink'
//     input:
//         path weighted_roc

//     output:
//         path "*.svg"

//     script:
//     def otu = "${(bam.getSimpleName() =~ /.+_(.+)/)[0][1]}"
//     def curve_str = weighted_roc{"--curve" + it +"="}.join(" ")
//         chr_name_match = params.chromosomes.collect{"chr" + it + "="}.join("\n")
//         chr_prefix = params.chromosomes.join("_")
//     """
//     "${params.rtg_dir_path}/rtg" rocplot -b ${ground_truth_vcf} -c ${test_vcf} -t "${ref.getSimpleName()}.sdf" -o "vcf_benchmark_${test_vcf.getSimpleName()} --"
//     ${params.rtg_dir_path}/rtg-tools-3.13/rtg rocplot \
//     --curve /hps/nobackup/goldman/ipoetzsch/boss_sequence/work/ce/5a858c7aaa6f2b498b63e97fe69da8/vcf_benchmark_control_33_hg002_chr21/weighted_roc.tsv.gz=control_33,cov=20x \
//     --curve /hps/nobackup/goldman/ipoetzsch/boss_sequence/work/b3/436a080c2865ad6211551fd1b18aa2/vcf_benchmark_boss_88_hg002_chr21/weighted_roc.tsv.gz=boss_88,cov=20x \
//     --curve /hps/nobackup/goldman/ipoetzsch/boss_sequence/work/f3/0f20a00ee105a15632b9d9b91fcffd/vcf_benchmark_control_88_hg002_chr21/weighted_roc.tsv.gz=control_88,cov=55x \
//     --curve /hps/nobackup/goldman/ipoetzsch/boss_sequence/work/76/d64893288437cb2e2a171bf0fa056e/vcf_benchmark_control_289_hg002_chr21/weighted_roc.tsv.gz=control_289,cov=180x \
//     --palette=blind_8 \
//     --svg=${title} \
//     --plain
//     """
// }


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


        

        // pile_csv = pileup_2(pile_full, otu.first())

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

// workflow BENCHMARK_VCF {
//     take:
//         ref
//         ground_truth_vcf
//         test_vcf

//     main:
//         benchmark_variants(ground_truth_vcf, test_vcf, ref)
    
// }

// workflow benchmark {
//     // Test_vcf
//     test_vcf = channel.fromPath( "${params.test_vcf_base}" )
//     test_vcf
//         .map { tuple( it.simpleName, it ) }
//         .groupTuple()
//         .set { test_vcf_collection }

//     channel.fromPath( "${params.test_vcf_base_idx}" )
//         .map { tuple( it.simpleName, it ) }
//         .combine( test_vcf_collection, by: 0 )
//         .transpose( by: 2 )
//         .map { _tmp, idx, vcf -> tuple( vcf, idx ) }
//         .set{test_vcf_tuple}

//     // Sequence log
//     ground_truth_vcf = channel.of(tuple("${params.ground_truth}", "${params.ground_truth_idx}"))
//     // Index reference
//     ref = channel.fromPath("${params.ref}")
//     benchmark_variants(ground_truth_vcf.first(), test_vcf_tuple, ref.first())
// }

