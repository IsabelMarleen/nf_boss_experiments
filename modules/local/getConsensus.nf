process getConsensus {
    debug true
    clusterOptions '--mem=4G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_vi}", mode: 'symlink'
    input:
        tuple path(vcf), path(vcf_idx)
        path ref
    output:
        path '*.fa', emit:consensus
    script:
    """
    module load bcftools
    bcftools consensus -f $ref $vcf -o consensus_${vcf.getBaseName(2)}.fa
    """
}