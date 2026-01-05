#!/usr/bin/env nextflow

// Run using 'nextflow run boss_analyse.nf -with-conda -resume'
// Depends on a nextflow.config file that defines the parameters

params.base_in_dir = "${params.base_out_dir}/${params.exp_name}"

process indexPaf{
    clusterOptions '--mem=8G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:10:00'
    publishDir "${params.output_dir_debug}/index", mode: 'symlink'
    conda "${params.conda_envs}/simulation.yaml"
    input:
        file ref
    output:
        path 'ref.mmi', emit: ref_idx
        path '*.fa', emit: ref_unzipped
    script:
        """
        minimap2 -d ref.mmi $ref
        gunzip $ref -c > ${ref.getSimpleName()}.fa
        """
}

process mapPaf {
    tag "${input.getSimpleName()}"
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    publishDir "${params.output_dir_debug}/map_paf", mode: 'symlink'
    conda "${params.conda_envs}/simulation.yaml"
    input:
        file input
        file index
    output:
        path '*.paf', emit: ref_path
    script:
    """
    minimap2 -x map-ont -t 32 --secondary=no -c $index $input -o ${input.getSimpleName()}.paf
    """

}

process separateTarget {
    tag "${paf.getSimpleName()}"
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:10:00'
    publishDir "${params.output_dir_debug}/sep_target", mode: 'symlink'
    conda "${params.conda_envs}/simulation.yaml"
    input:
        tuple file(paf), file(reads)
    output:
        path '*.fq', emit: fq
    script:
    """
    python3 ${params.analyse_script_path}/separate_by_target.py $paf $reads
    """

}

process mapSam {
    tag "${reads.getSimpleName()}"
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    publishDir "${params.output_dir_debug}/map_sam", mode: 'symlink'
    conda "${params.conda_envs}/simulation.yaml"
    input:
        file reads
        file unzipped_ref
    output:
        tuple(path('*.bam'), path('*.bam.bai'), emit: otu_mapped)
    script:
    def otu = "${(reads.getSimpleName() =~ /.+_(.+)/)[0][1]}"
    """
    python3 ${params.analyse_script_path}/extract_sequence.py ${unzipped_ref} "${otu}" > ${params.exp_name}_${otu}.fa &&
    minimap2 -ax map-ont -t 32 --secondary=no --sam-hit-only -c ${params.exp_name}_${otu}.fa ${reads} |
    samtools sort -@ 32 | samtools view -b > ${reads.getSimpleName()}.bam &&
    samtools index -@ 32 ${reads.getSimpleName()}.bam
    """
}

process pileup {
    tag "${bam.getSimpleName()}"
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    publishDir "${params.output_dir_debug}/pileup", mode: 'symlink'
    conda "${params.conda_envs}/simulation.yaml"
    input:
        tuple path(bam), path(bai)
    output:
        path '*.pup', emit: pup
        path "${bam.getSimpleName()}.csv", emit:csv
        path '*_full.csv'
    script:
    def otu = "${(bam.getSimpleName() =~ /.+_(.+)/)[0][1]}"
    """
    samtools mpileup -Q 0 ${bam} > ${bam.getSimpleName()}.pup &&
    python3 ${params.analyse_script_path}/process_pileup.py ${bam.getSimpleName()}.pup ${otu} ${bam.getSimpleName()}_full.csv > ${bam.getSimpleName()}.csv
    """
}

process recordUnblocks {
    tag "${fq.getSimpleName()}"
    clusterOptions '--mem=2G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_debug}/unblocks", mode: 'symlink'
    conda "${params.conda_envs}/simulation.yaml"
    input:
        path fq
    output:
        path '*.csv', emit: csv
    script:
    """
    python3 ${params.analyse_script_path}/record_unblocks.py $fq ${params.mu} > ${fq.getSimpleName()}_unblocks.csv
    """
}


process createUnblock_dataframe {
    clusterOptions '--mem=4G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_debug}/results", mode: 'symlink'
    input:
        path csv
    output:
        path 'unblocks.csv', emit:unblocks
    script:
    """
    cat <(echo cond,time,otu,total,base_total,unb,unb_ratio) $csv > unblocks.csv
    """
}

process create_coverage_dataframe { 
    clusterOptions '--mem=4G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_debug}/results", mode: 'symlink'
    input:
        path csv
    output:
        path 'coverage.csv', emit: coverage

    script: 
    """
    mkdir -p results/$params.exp_name &&
    cat <(echo cond,time,otu,mean_coverage,low_coverage_prop,evenness) $csv > coverage.csv
    """
} 

// Analyse read log
process analyse_log {
    clusterOptions '--mem=8G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_debug}/results", mode: 'symlink'
    conda "${params.conda_envs}/simulation.yaml"
    input: 
        path log_file
    output: 
        path "*.csv", emit: analysed_log

    script:
    """
    python3 ${params.analyse_script_path}/analyse_boss_log.py $log_file ${log_file.getSimpleName()}_log_processed.csv
    """
}

// Visualise results
process visualise_simulation {
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:20:00'
    publishDir "${params.output_dir_debug}/results", mode: 'copy'
    conda "${params.conda_envs}/visualisation.yaml"
    cache false
    input: 
        path coverage
        path unblocks
        path analysed_log
    output: 
        path "*.pdf", emit: result_plot

    script:
    """
    #now=date +'%Y/%m/%d'
    mkdir -p ${params.output_dir_results} &&
    Rscript ${params.analyse_script_path}/visualise_simulation.R \
    --input_cov ${coverage} \
    --input_unb ${unblocks} \
    --dump_time ${params.dump_time} \
    --pores ${params.pores} \
    --analysed_log $analysed_log \
    --output ${params.exp_name}_plots.pdf
    """
}

process benchmark_variants {
    tag "${test_vcf.getSimpleName()}"
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:20:00'
    publishDir "${params.output_dir_debug}/results", mode: 'symlink'
    input:
        tuple path(ground_truth_vcf), path(ground_truth_vcf_idx)
        tuple path(test_vcf), path(test_vcf_idx)
        path ref

    output:
        path "vcf_benchmark_*"
        path "*.sdf"

    script:
    """
    "${params.rtg_dir_path}/rtg" format -o "${ref.getSimpleName()}.sdf" ${ref}
    "${params.rtg_dir_path}/rtg" vcfeval -b ${ground_truth_vcf} -c ${test_vcf} -t "${ref.getSimpleName()}.sdf" -o "vcf_benchmark_${test_vcf.getSimpleName()}"
    """
}

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
            .map { tuple( it.simpleName, it ) }
            .groupTuple()
            .set { input_reads_collection }

        // Index reference
        processed_ref = indexPaf(ref_input)
        ref_idx = processed_ref.ref_idx
        ref_unzipped = processed_ref.ref_unzipped

        // Map paf
        first_map = mapPaf(input_reads, ref_idx.first())
                
        first_map
                .map { tuple( it.simpleName, it ) }
                .combine( input_reads_collection, by: 0 )
                .transpose( by: 2 )
                .map { _case_id, paf, fa -> tuple( paf, fa ) }
                .set{first_map_tuple}
                

        // Separate target
        sep_target = separateTarget(first_map_tuple)

        // Map sam
        mapped_sam = mapSam(sep_target, ref_unzipped.first())

        // Pileup
        pile = pileup(mapped_sam)

        // Record unblocks
        unblocks_ind = recordUnblocks(sep_target)

        // Create unblocks dataframe
        unblocks = createUnblock_dataframe(unblocks_ind.collect())

        // Create coverage dataframe
        cov = create_coverage_dataframe(pile.csv.collect())

        // Analyse log
        analysed_log = analyse_log(seq_log)

        // // Create plots
        plots = visualise_simulation(cov, unblocks, analysed_log)
    emit:
        coverage = cov
        unblocks = unblocks
        plots = plots
}

workflow BENCHMARK_VCF {
    take:
        ref
        ground_truth_vcf
        test_vcf

    main:
        benchmark_variants(ground_truth_vcf, test_vcf, ref)
    
}

workflow benchmark {
    // Test_vcf
    test_vcf = channel.fromPath( "${params.test_vcf_base}" )
    test_vcf
        .map { tuple( it.simpleName, it ) }
        .groupTuple()
        .set { test_vcf_collection }

    channel.fromPath( "${params.test_vcf_base_idx}" )
        .map { tuple( it.simpleName, it ) }
        .combine( test_vcf_collection, by: 0 )
        .transpose( by: 2 )
        .map { _tmp, idx, vcf -> tuple( vcf, idx ) }
        .set{test_vcf_tuple}

    // Sequence log
    ground_truth_vcf = channel.of(tuple("${params.ground_truth}", "${params.ground_truth_idx}"))
    // Index reference
    ref = channel.fromPath("${params.ref}")
    vars = benchmark_variants(ground_truth_vcf.first(), test_vcf_tuple, ref.first())
    
    emit:
        variants = vars
}

workflow  {
    // Input reads
    input_reads = channel.fromPath( "${params.reads}*.fa" )
    // Sequence log
    seq_log = channel.of("${params.log}")
    // Index reference
    ref_input = channel.fromPath("${params.ref}")

    ANALYSE_BOSS(input_reads, seq_log, ref_input)
}