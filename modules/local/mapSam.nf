process mapSam {
    tag "${reads.getSimpleName()}"
    clusterOptions '--mem=32G --nodes=1 --cpus-per-task=32 --ntasks=1 --time=00:30:00'
    conda "${params.conda_envs}/simulation.yaml"
    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/90/90b05084d03e6b3d92e923b983f0ed5785b03169521d4e3bc85152f6c3b32122/data"
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