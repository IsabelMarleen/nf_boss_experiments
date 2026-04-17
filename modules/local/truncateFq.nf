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