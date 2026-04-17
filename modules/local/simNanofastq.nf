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
    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7f795a3f3858f4e53d9c2d615e66f704bf27ad6cf1b291901694143a1bd31294/data"
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
    simulator.py genome -rg ${consensus_ref} -c ${params.model} -n $readnumber --fastq -t 32 -min 401 --seed 4444600 -o ${consensus_ref.getBaseName()}
    """
} 