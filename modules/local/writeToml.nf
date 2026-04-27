// Create toml file from preprocessing
process writeToml {
    input: 
        path full_fq
        path full_paf
        path trunc_paf
        path trunc_fq
        path fq_offsets
        path full_paf_offsets
        path trunc_paf_offsets
        path reference
    output: 
        path("*.toml")

    script:
    maxb = params.readnumber.values().sum()*0.6/params.batch_size as Integer
    if (params.barcodes == null){
        bc_string = ""
    }
    else{
        bc = params.barcodes.values().join('\\", \\"')
        bc_string = "\nbarcodes = "+'[\\"' + bc + '\\"]'
    }
    """
    echo "[general]
name = '"${params.exp_name}"'                   # experiment name
ref = '"\$PWD/${reference}"'                        # reference fasta file. Not specifying a file switches to BOSS-AEONS
wait = 30                       # waiting time between updates in live version

[simulation]
fq = '"\$PWD/${full_fq}"'                   # fastq file
paf_full = '"\$PWD/${full_paf}"'     # full mapping paf file
paf_trunc = '"\$PWD/${trunc_paf}"'            # truncated mapping paf file
maxb = $maxb                  # maximum number of batches
batchsize = $params.batch_size # How many reads are processed at once
binit = 10
dumptime = ${params.dump_time} # how much pseudotimes needs to be surpassed to write to file -- adjusted to better match 35 Mio. from before pseudotime adjustment
accept_unmapped = ${params.accept_unmapped}${bc_string}" > "${params.exp_name}.toml"
    """
}