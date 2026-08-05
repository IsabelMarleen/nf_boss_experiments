process getGunzipReads {
    input:
        val link
    output:
        path '*.fastq'
    script:
    """
    wget ${link}
    gunzip -k *.gz
    """

}