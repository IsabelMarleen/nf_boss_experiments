process simNanofastq { 
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