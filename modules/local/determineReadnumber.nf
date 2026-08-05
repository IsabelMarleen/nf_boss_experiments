process determineReadnumber {
    input:
        path fastq
    output:
        env 'float'
    script:
    """
    float=`wc -l ${fastq} | awk '{print \$1 / 4 }'`
    """
}