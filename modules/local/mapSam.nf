process mapSam {
    tag "${reads.getSimpleName()}"
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    conda "${params.conda_envs}/simulation.yaml"
    input:
        file reads
        file unzipped_ref
        file otu_file
    output:
        tuple(path('*.bam'), path('*.bam.bai'), emit: otu_mapped)
    script:
    def otu = "${(reads.getSimpleName() =~ /.+_(.+)/)[0][1]}"
    """
    ${params.script_dir}extract_sequence.py ${unzipped_ref} "${otu}" > ${params.exp_name}_${otu}.fa &&
    minimap2 -ax map-ont -t 32 --secondary=no --sam-hit-only -c ${params.exp_name}_${otu}.fa ${reads} |
    samtools sort -@ 32 | samtools view -b > ${reads.getSimpleName()}.bam &&
    samtools index -@ 32 ${reads.getSimpleName()}.bam
    """
}