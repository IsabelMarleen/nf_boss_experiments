#!/usr/bin/env nextflow

// Run using 'nextflow run boss_analyse.nf -with-conda -resume -c ../nextflow.config'

params.exp_name = "benchmark_hg002_chr21"

params.base_in_dir = "${params.base_out_dir}/${params.exp_name}"
// params.br_output = "${params.base_out_preprocess}/br_output"


params.log_file = "20250905-193019_boss.log"

// params.path = "${params.br_output}/boss_runs_out"
// params.reads = "${params.path}/${params.exp_name}/00_reads/"
params.reads = "/hps/nobackup/goldman/ipoetzsch/boss_experiments/work/97/7fed3d495da0b6b53cde1208eb9d9c/00_reads/"
params.log = "/hps/nobackup/goldman/ipoetzsch/boss_experiments/work/97/7fed3d495da0b6b53cde1208eb9d9c/20251123-133405_boss.log"
// params.log = "${params.path}/${params.exp_name}/${params.log_file}"
// params.ref = "${params.base_out_dir}/${params.exp_name}/variant_input/Homo_sapiens.GRCh38.dna.chromosome.21.fa"
params.ref = "/hps/nobackup/goldman/ipoetzsch/boss_experiments/work/a0/c4632f88faa413f3b321fe0c40c1ce/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz"
params.analyse_script_path = "${projectDir}/scripts"
params.conda = "${params.conda_base_dir}"

params.mu = 400
params.dump_time = 35_000_000
params.pores = 512

params.rtg_dir_path = "${params.software_dir}/rtg-tools-3.13/"

params.base_out_dir = "${params.base_out_dir}/${params.exp_name}"
params.output_dir_results = "${params.base_out_dir}/results"
params.output_dir_debug = "${params.base_out_dir}/debug"


process indexPaf{
    clusterOptions '--mem=8G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:10:00'
    publishDir "${params.output_dir_debug}/index", mode: 'symlink'
    conda "${params.conda}/boss_simulation"
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
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    publishDir "${params.output_dir_debug}/map_paf", mode: 'symlink'
    conda "${params.conda}/boss_simulation"
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
    clusterOptions '--mem=2G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_debug}/sep_target", mode: 'symlink'
    conda "${params.conda}/boss_simulation"
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
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    publishDir "${params.output_dir_debug}/map_sam", mode: 'symlink'
    conda "${params.conda}/boss_simulation"
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
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    publishDir "${params.output_dir_debug}/pileup", mode: 'symlink'
    conda "${params.conda}/boss_simulation"
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
    clusterOptions '--mem=2G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_debug}/unblocks", mode: 'symlink'
    conda "${params.conda}/boss_simulation"
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
    conda "${params.conda}/boss_simulation"
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
    conda "${params.conda}/boss_visualise"
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

// // Convert to bam
// process getBam {
//     clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=03:00:00'
//     publishDir "${params.output_dir_br_output}", mode: 'symlink'
//     conda "${params.conda}/boss4"
//     input: 
//         path fastq
//     output: 
//         path "*.bam", emit: bam

//     script:
//     """
//     module load samtools
//     minimap2 -x -a map-ont -t 32 --secondary=no -c $params.ref $fastq | \
//     samtools fixmate -u -m - - | \
//     samtools sort -u -@2 -T /tmp/example_prefix - | \
//     samtools markdup -@8 --reference $params.ref -o ${fastq.getSimpleName()}.bam

//     """
// }

process benchmark_variants {
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:20:00'
    publishDir "${params.output_dir_debug}/results", mode: 'symlink'
    input:
        path ground_truth_vcf
        path test_vcf
        path ref

    output:
        path "vcf_benchmark_*"
        path "*.sdf"

    script:
    """
    ".${params.rtg_dir_path}/rtg" format -o "${ref.getSimplename()}.sdf" ${ref}
    ".${params.rtg_dir_path}/rtg" vcfeval -b ${ground_truth_vcf} -c ${test_vcf} -t "${ref.getSimplename()}.sdf" -o "vcf_benchmark_${test_vcf.getsimpleName()}"
    """
}

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
                .map { case_id, paf, fa -> tuple( paf, fa ) }
                .set{first_map_tuple}
                

        // Separate target
        sep_target = separateTarget(first_map_tuple)

        // Map sam
        mapped_sam = mapSam(sep_target, ref_unzipped)

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

workflow  {
    // Input reads
    input_reads = channel.fromPath( "${params.reads}*.fa" )
    // Sequence log
    seq_log = channel.of("${params.log}")
    // Index reference
    ref_input = channel.fromPath("${params.ref}")

    ANALYSE_BOSS(input_reads, seq_log, ref_input)
}