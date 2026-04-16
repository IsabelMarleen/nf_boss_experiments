process pileup {
    debug true
    tag "${meta.join("_")}"
    clusterOptions '--nodes=1 --cpus-per-task=32 --ntasks=1'
    memory { 1.5.GB * bam.size() * task.attempt }
    time {1.min * bam.size() * task.attempt}
    maxRetries 3
    conda "${params.conda_envs}/simulation.yaml"
    errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
    input:
        tuple val(meta), path(bam)
        tuple val(meta2), path(bai)
        file otu_file
    output: 
        // path '*.pup', emit: pup
        // path "${bam.getSimpleName()}.csv", emit:csv
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
    ${params.script_dir}process_pileup.py ${meta.join("_")}.pup $otu ${meta.join("_")}.csv ${bam.size()}
    """
}