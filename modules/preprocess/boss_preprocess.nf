#!/usr/bin/env nextflow

// Run using 'nextflow run boss_preprocess.nf -entry PREPROCESS_BOSS -with-conda -resume -c ../nextflow.config'

/*
 * Pipeline parameters
 */
params.exp_name = "benchmark_hg002_chr21"

// params.input_base_dir -- defined in nextflow.config
params.genome = "${params.base_input_dir}/${params.exp_name}/variant_input/consensus_21_HG002_GRCh38_1_22_v4.2.1_benchmark.fa"
params.ref = "${params.base_input_dir}/${params.exp_name}/variant_input/Homo_sapiens.GRCh38.dna.chromosome.21.fa"

params.abs_path_to_nanosim_repo = "${params.software_dir}/NanoSim"
params.model = "${params.abs_path_to_nanosim_repo}/pre-trained_models/human_giab_hg002_sub1M_kitv14_dorado_v3.2.1/training"
params.readnumber = 1000000

params.abspath_to_boss_runs_repo = "${params.software_dir}/BOSS-RUNS2"
params.mu = 400
params.dump_time = 35000000
params.accept_unmapped = "true"

params.base_out_preprocess = "${params.base_out_dir}/${params.exp_name}"
params.output_dir_nanosim = "${params.base_out_preprocess}/nanosim_out"
params.output_dir_br_input = "${params.base_out_preprocess}/br_input"
params.output_dir_br_output = "${params.base_out_preprocess}/br_output"


process simNanofastq { 
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    publishDir (
        path: "${params.output_dir_nanosim}", 
        mode: 'symlink',
        saveAs: {fn ->
            if (fn.endsWith(".fastq")) { "${file(fn).getBaseName(1)}.fq" }
            else {"${fn}" }
        }
        
    )
    conda "${params.conda_base_dir}/nanosim"
    input:
        path consensus_ref

    output:
        path '*_aligned_reads.fastq', emit: aligned_fastq
        path '*_unaligned_reads.fastq', emit: fastq
        path '*error_profile', emit: error_profile

    script: 
    """
    $params.abs_path_to_nanosim_repo/src/simulator.py genome -rg ${consensus_ref} -c ${params.model} -n ${params.readnumber} --fastq -t 32 -min 401 --seed 4444600

    """
} 

// Truncate fq Step 1 of prepare input for BRsim
process truncateFq {
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:20:00'
    publishDir(
        path: "${params.output_dir_br_input}", 
        mode: 'symlink',        
        saveAs: {fn ->
            if (fn.endsWith(".fastq")) { "${file(fn).getBaseName(1)}.fq" }
            else {"${fn}" }
        }
    )
    conda "${params.conda_base_dir}/boss2"
    input: 
        path aligned_fastq
    output: 
        path "*_trunc.fq", emit:trunc_fq
        path aligned_fastq, emit:full_fq

    script:
    """
    cut -c -${params.mu} ${aligned_fastq} > ${aligned_fastq.getSimpleName()}_trunc.fq
    """
}

// Scan fq offsets Step 2 of prepare input for BRsim
process scan_offsetsFq {
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    publishDir "${params.output_dir_br_input}", mode: 'symlink'
    conda "${params.conda_base_dir}/boss2"
    input: 
        path aligned_fastq
    output: 
        path "*.offsets.npy", emit: fq_offsets

    script:
    """
    python3 ../../scripts/scan_offsets_fq_nf.py $aligned_fastq
    """
}

// Map paf Step 3 of prepare input for BRsim
process mapPaf {
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    publishDir "${params.output_dir_br_input}", mode: 'symlink'
    conda "${params.conda_base_dir}/boss2"
    input: 
        path aligned_fastq
        path ref
    output: 
        path "*.paf", emit: mappings

    script:
    """
    # mm2 = which("minimap2", path='/'.join(executable.split('/')[0:-1])).strip()
    minimap2 -x map-ont -t 32 --secondary=no -c $ref $aligned_fastq > ${aligned_fastq.getSimpleName()}.paf
    """
}

// Map paf trunc Step 4 of prepare input for BRsim
process mapPaf_trunc {
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    publishDir "${params.output_dir_br_input}", mode: 'symlink'
    conda "${params.conda_base_dir}/boss2"
    input: 
        path aligned_fastq_trunc
        path ref
    output: 
        path "*.paf", emit: mappings_trunc

    script:
    """
    # from shutil import which

    # mm2 = which("minimap2", path='/'.join(executable.split('/')[0:-1])).strip()
    minimap2 -x map-ont -t 32 --secondary=no -c $ref $aligned_fastq_trunc > ${aligned_fastq_trunc.getSimpleName()}.paf
    """
}

// Scan paf offsets Step 5 of prepare input for BRsim
process scan_offsets_Paf {
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    publishDir "${params.output_dir_br_input}", mode: 'symlink'
    conda "${params.conda_base_dir}/boss2"
    input: 
        path mappings
        path mappings_trunc
    output: 
        path "${mappings.getSimpleName()}.paf.offsets", emit: mappings_offsets
        path "*_trunc.paf.offsets", emit: mappings_offsets_trunc

    script:
    """
    python3 ../../scripts/scan_offsets_paf_nf.py ${mappings} ${mappings_trunc}
    """
}

// Create toml file from preprocessing
process writeToml {
    clusterOptions '--mem=1G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_br_input}", mode: 'copy'
    // cache false
    input: 
        path full_fq
        path full_paf
        path trunc_paf
        path trunc_fq
        path fq_offsets
        path full_paf_offsets
        path trunc_paf_offsets
        path reference
    output: 
        path("*.toml")

    script:
    maxb = params.readnumber*0.6/1000 as Integer
    """
    echo "[general]
name = '"${params.exp_name}"'                   # experiment name
ref = '"${reference}"'                        # reference fasta file. Not specifying a file switches to BOSS-AEONS

[simulation]
fq = '"${full_fq}"'                   # fastq file
paf_full = '"${full_paf}"'     # full mapping paf file
paf_trunc = '"${trunc_paf}"'            # truncated mapping paf file
maxb = $maxb                  # maximum number of batches
batchsize = 1000 # How many reads are processed at once
binit = 10
dumptime = ${params.dump_time} # how much pseudotimes needs to be surpassed to write to file -- adjusted to better match 35 Mio. from before pseudotime adjustment
accept_unmapped = ${params.accept_unmapped}" > "${params.exp_name}.toml"
    """
}

// Run br sim
process runBRSim {
    clusterOptions '--mem=32G --nodes=4 --cpus-per-task=32 --ntasks=1 --time=03:00:00'
    publishDir "${params.output_dir_br_output}", mode: 'symlink'
    conda "${params.conda_base_dir}/boss4"
    input: 
        path toml
    output: 
        path("out_*"), emit: dir
        path("00_reads"), emit: reads_dir
        path("*_boss.log"), emit: log

    script:
    """
    boss --toml ${params.toml}
    """
}

// // Convert to bam
// process getBam {
//     clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=03:00:00'
//     publishDir "${params.output_dir_br_output}", mode: 'symlink'
//     conda "${params.conda_base_dir}/boss4"
//     input: 
//         path fastq
//         path toml
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

workflow PREPROCESS_BOSS{
    take:
        input_genome
        ref
    main:
        // Simulate reads using Nanosim
        nanosim = simNanofastq(input_genome)

        // Prepare input for BRsim
        // Truncate fq Step 1 of prepare input for BRsim
        trunc_fq = truncateFq(nanosim.aligned_fastq)

        // Scan fq offsets Step 2 of prepare input for BRsim
        fq_combined = trunc_fq.full_fq.mix(trunc_fq.trunc_fq)
        fq_offsets = scan_offsetsFq(fq_combined)

        // Map paf Step 3 of prepare input for BRsim
        mpaf = mapPaf(nanosim.aligned_fastq, input_genome)

        // Map paf trunc Step 4 of prepare input for BRsim
        mpaf_trunc = mapPaf_trunc(trunc_fq.trunc_fq, input_genome)

        // Scan paf offsets Step 5 of prepare input for BRsim
        paf_offsets = scan_offsets_Paf(mpaf.mappings, mpaf_trunc.mappings_trunc)

        // Create toml file based on input files above
        toml = writeToml(trunc_fq.full_fq, mpaf, mpaf_trunc, trunc_fq.trunc_fq, fq_offsets.fq_offsets.collect(), 
                        paf_offsets.mappings_offsets, paf_offsets.mappings_offsets_trunc, ref)


    emit:
        toml
}

workflow  PREPROCESS_WITH_SEQ{
    input_genome = Channel.of(params.genome)
    reference = Channel.of(params.ref)
    // Run preprocessing
    toml = PREPROCESS_BOSS(input_genome, reference) 

    // Run brsim
    runBRSim(toml)
}

workflow  {
    // Download reference file(s) for each chromosome and benchmark vcfs + index 
    input_genome = Channel.of(params.genome)
    reference = Channel.of(params.ref)

    PREPROCESS_BOSS(input_genome, reference)
}