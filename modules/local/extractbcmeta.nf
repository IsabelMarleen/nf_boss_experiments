process extract_bc_meta{
    tag "${file.getSimpleName()}"
    input:
        path file
    output:
        tuple val(meta), path(file)
    script:
        labels = params.barcodes.keySet().join("|")
        def matcher = file =~ "$labels"
        long_bc = params.barcodes[matcher[0]]
        subs = long_bc.substring(7) as Integer as String
        meta = [bc: "bc".plus(subs)]
    """
    touch ${file}
    """
}