process mapSam {
    tag "${reads.getSimpleName()}"
    input:
        file reads
        file unzipped_ref
        file otu_file
    output:
        tuple(path('*.bam'), path('*.bam.bai'), emit: otu_mapped)
    script:
    def otu = "${(reads.getSimpleName() =~ /.+_(.+)/)[0][1]}"
    """
    extract_sequence.py ${unzipped_ref} "${otu}" > ${params.exp_name}_${otu}.fa &&
    minimap2 -ax map-ont -t 32 --secondary=no --sam-hit-only -c ${params.exp_name}_${otu}.fa ${reads} |
    samtools sort -@ 32 | samtools view -b > ${reads.getSimpleName()}.bam &&
    samtools index -@ 32 ${reads.getSimpleName()}.bam
    """
}