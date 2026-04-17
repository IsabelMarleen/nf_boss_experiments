process pileup {
    debug true
    tag "${meta.join("_")}"
    clusterOptions '--nodes=1 --cpus-per-task=32 --ntasks=1'
    memory { 1.5.GB * bam.size() * task.attempt }
    time {1.min * bam.size() * task.attempt}
    maxRetries 3
    conda "${params.conda_envs}/simulation.yaml"
    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b6/b62d54a49751d6761d41b48f7b81b8450190946068bfef21d83a165fd48e4c54/data"
    errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
    input:
        tuple val(meta), path(bam)
        tuple val(meta2), path(bai)
        file otu_file
    output: 
        path 'tmp.csv', emit: pup
        path "${meta.join("_")}.csv", emit: csv
    script:
    if (params.barcodes == null){
        otu = meta.get(1)
        assert otu == meta2.get(1)
    }
    else{
        otu = meta.get(2)
        assert otu == meta2.get(2)
    }
    """
    samtools mpileup -Q 0 ${bam} > ${meta.join("_")}.pup &&
    process_pileup.py ${meta.join("_")}.pup $otu ${meta.join("_")}.csv ${bam.size()}
    """
}