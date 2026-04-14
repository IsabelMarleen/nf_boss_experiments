#!/usr/bin/env nextflow

// Run using 'nextflow run boss_preprocess.nf -entry PREPROCESS_BOSS -with-conda -resume -c ../nextflow.config'

/*
 * Pipeline parameters in ../nextflow.config
 */

params.script_path = "/hps/software/users/goldman/ipoetzsch/nf_boss_experiments/bin"
process simNanofastq { 
    clusterOptions '--nodes=1 --cpus-per-task=32 --ntasks=1' //--mem=32G  --time=00:30:00
    memory { 32.GB * params.chromosomes.size() * task.attempt }
    time { 1.h * params.chromosomes.size() * task.attempt }
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
    labels = params.readnumber.keySet().join("|")

    def matcher = consensus_ref =~ "$labels"
    readnumber = params.readnumber[matcher[0]] as Integer
    """
    $params.abs_path_to_nanosim_repo/src/simulator.py genome -rg ${consensus_ref} -c ${params.model} -n $readnumber --fastq -t 32 -min 401 --seed 4444600 -o ${consensus_ref.getBaseName()}
    """
} 

process addBarcodes{
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=8 --ntasks=1 --time=00:20:00'
    publishDir(
        path: "${params.output_dir_br_input}", 
        mode: 'symlink'
    )
    input: 
        path simulated_fastq
    output: 
        path "*_barcoded.fq", emit:barcoded_fq

    script:
    labels = params.barcodes.keySet().join("|")
    def matcher = simulated_fastq =~ "$labels"
    barcode = params.barcodes[matcher[0]]
    """
    sed '1~4 s/\$/ barcode=$barcode/' ${simulated_fastq} > ${simulated_fastq.getSimpleName()}_barcoded.fq
    """
}

process combineBarcodedFq{
    clusterOptions '--nodes=1 --cpus-per-task=8 --ntasks=1'
    memory { 64.GB * params.readnumber.values().size() * task.attempt }
    time { 20.min * params.readnumber.values().size() * task.attempt }
    maxRetries 3
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }
    publishDir(
        path: "${params.output_dir_br_input}", 
        mode: 'symlink'
    )
    input: 
        path barcoded_fastq
    output: 
        path "*_barcoded.fq", emit:combined_barcoded_fq

    script:
    """
    cat ${barcoded_fastq} | \
    awk '{OFS="\\t"; getline seq; \
                getline sep; \
                getline qual; \
                print \$0,seq,sep,qual}' | \
    shuf | \
    awk '{OFS="\\n"; print \$1 " " \$2,\$3,\$4,\$5}' \
    > ${params.exp_name}_combined_barcoded.fq
    """
}

// Truncate fq Step 1 of prepare input for BRsim
process truncateFq {
    clusterOptions '--nodes=1 --cpus-per-task=8 --ntasks=1'
    memory { 32.GB * params.readnumber.values().size() * task.attempt }
    time { 40.min * params.readnumber.values().size() * task.attempt }
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
        path sampler_func_file, stageAs: 'sampler.py'

    output: 
        path "*.offsets.npy", emit: fq_offsets

    script:
    """
    scan_offsets_fq_nf.py $aligned_fastq
    """
}

// Map paf Step 3 of prepare input for BRsim
process mapPaf {
    clusterOptions '--nodes=1 --cpus-per-task=32 --ntasks=1'
    memory { 32.GB * params.readnumber.values().size() * params.chromosomes.size() * task.attempt }
    time {20.min * params.readnumber.values().size() * params.chromosomes.size() * task.attempt}
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
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1'
    time {30.min * params.readnumber.values().size() * params.chromosomes.size() * task.attempt}
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
    scan_offsets_paf_nf.py ${mappings} ${mappings_trunc}
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
    maxb = params.readnumber.values().sum()*0.6/params.batch_size as Integer
    if (params.barcodes == null){
        bc_string = ""
    }
    else{
        bc = params.barcodes.values().join('\\", \\"')
        bc_string = "\nbarcodes = "+'[\\"' + bc + '\\"]'
    }
    """
    echo "[general]
name = '"${params.exp_name}"'                   # experiment name
ref = '"\$PWD/${reference}"'                        # reference fasta file. Not specifying a file switches to BOSS-AEONS
wait = 30                       # waiting time between updates in live version

[simulation]
fq = '"\$PWD/${full_fq}"'                   # fastq file
paf_full = '"\$PWD/${full_paf}"'     # full mapping paf file
paf_trunc = '"\$PWD/${trunc_paf}"'            # truncated mapping paf file
maxb = $maxb                  # maximum number of batches
batchsize = $params.batch_size # How many reads are processed at once
binit = 10
dumptime = ${params.dump_time} # how much pseudotimes needs to be surpassed to write to file -- adjusted to better match 35 Mio. from before pseudotime adjustment
accept_unmapped = ${params.accept_unmapped}${bc_string}" > "${params.exp_name}.toml"
    """
}

// Run br sim
process runBRSim {
    clusterOptions '--mem=128G --nodes=4 --cpus-per-task=64 --ntasks=1 --time=12:00:00'
    publishDir "${params.output_dir_br_output}/11-22_20-53_rerun", mode: 'symlink'
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

workflow PREPROCESS_BOSS{
    take:
        input_genome
        ref
    main:
        // Simulate reads using Nanosim
        nanosim = simNanofastq(input_genome)

        // If barcoded experiment, add barcode info to each fq and combine
        if (params.barcodes != null){
            barcoded = addBarcodes(nanosim.aligned_fastq)
            aligned_fastq = combineBarcodedFq(barcoded.collect())
        }
        else{
            aligned_fastq = nanosim.aligned_fastq
        }

        // Prepare input for BRsim
        // Truncate fq Step 1 of prepare input for BRsim
        trunc_fq = truncateFq(aligned_fastq)

        // Scan fq offsets Step 2 of prepare input for BRsim
        boss_sampler = Channel.fromPath(params.abs_path_to_boss_repo+"/sampler.py")
        fq_combined = trunc_fq.full_fq.mix(trunc_fq.trunc_fq)
        fq_offsets = scan_offsetsFq(fq_combined, boss_sampler)
        // TODO: Combine input genomes into one file using cat

        // Map paf Step 3 of prepare input for BRsim
        mpaf = mapPaf(aligned_fastq, ref)

        // Map paf trunc Step 4 of prepare input for BRsim
        mpaf_trunc = mapPaf_trunc(trunc_fq.trunc_fq, ref)

        // Scan paf offsets Step 5 of prepare input for BRsim
        paf_offsets = scan_offsets_Paf(mpaf.mappings, mpaf_trunc.mappings_trunc)

        // Create toml file based on input files above
        toml = writeToml(trunc_fq.full_fq, mpaf, mpaf_trunc, trunc_fq.trunc_fq, fq_offsets.fq_offsets.collect(), 
                        paf_offsets.mappings_offsets, paf_offsets.mappings_offsets_trunc, ref)

        // runBRSim(toml)
    emit:
        toml
}


workflow PREPROCESS_From_seq{
    take:
        full_fq
        ref
    main:
        // Prepare input for BRsim
        // Truncate fq Step 1 of prepare input for BRsim
        trunc_fq = truncateFq(full_fq)

        // Scan fq offsets Step 2 of prepare input for BRsim
        boss_sampler = Channel.fromPath(params.abs_path_to_boss_repo+"/sampler.py")
        fq_combined = full_fq.mix(trunc_fq.trunc_fq)
        fq_offsets = scan_offsetsFq(fq_combined, boss_sampler)

        // Map paf Step 3 of prepare input for BRsim
        mpaf = mapPaf(full_fq, ref)

        // Map paf trunc Step 4 of prepare input for BRsim
        mpaf_trunc = mapPaf_trunc(trunc_fq.trunc_fq, ref)

        // Scan paf offsets Step 5 of prepare input for BRsim
        paf_offsets = scan_offsets_Paf(mpaf.mappings, mpaf_trunc.mappings_trunc)

        // Create toml file based on input files above
        toml = writeToml(full_fq, mpaf, mpaf_trunc, trunc_fq.trunc_fq, fq_offsets.fq_offsets.collect(), 
                        paf_offsets.mappings_offsets, paf_offsets.mappings_offsets_trunc, ref)

        // runBRSim(toml)
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

// workflow  {
//     // Download reference file(s) for each chromosome and benchmark vcfs + index 
//     input_genome = Channel.of(params.genome)
//     reference = Channel.of(params.ref)

//     PREPROCESS_BOSS(input_genome, reference)
// }

workflow  {
    // Download reference file(s) for each chromosome and benchmark vcfs + index 
    reference = Channel.fromPath("${params.abs_path_to_boss_repo}/data/zymo.fa")
    full_fq = Channel.fromPath("${params.software_dir}/boss_exp/barcoded_input_data/FAT91932_pass_e7bf7751_f43c451e_4_bc.fastq")

    PREPROCESS_From_seq(full_fq, reference)
}