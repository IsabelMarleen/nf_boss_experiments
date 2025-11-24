#!/usr/bin/env nextflow

// Run using 'nextflow run boss_variant_input.nf -entry VARIANT_INPUT -with-conda -resume -c ../nextflow.config'

/*
 * Pipeline parameters
 */
params.exp_name = "benchmark_hg002_chr21"
// params.chromosomes = ["21", "22"]
params.chromosomes = "21"
// params.link_ref_base = "https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome."//${params.chr}.fa.gz"
params.link_ref = "https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome."//${params.chromosomes}.fa.gz"

params.link_vcf = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
params.link_vcf_idx = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"

// params.base_out_dir is defined in the config file
params.output_dir_vi = "${params.base_out_dir}/${params.exp_name}/variant_input"


process getRef {
    clusterOptions '--mem=2G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_vi}", mode: 'symlink'
    input:
        val ref_link
        val chr
    output:
        path '*.fa.gz', emit: ref_path
    script:
    """
    wget "${ref_link}${chr}.fa.gz"
    """

}

process combineRef{
    clusterOptions '--mem=2G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_vi}", mode: 'symlink'
    input:
        path refs
    output:
        path '*.fa.gz', emit: ref_path
    script:
    """
    # cat "${refs}" > "Homo_sapiens.GRCh38.dna.chromosome.${params.chromosomes.join(".")}.fa.gz"
    cat "${refs}" > "Homo_sapiens.GRCh38.dna.chromosome.${params.chromosomes}.fa.gz"
    """
}

process getVCF {
    clusterOptions '--mem=2G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_vi}", mode: 'symlink'
    input:
        val vcf_link
        val idx_link
    output:
        tuple path('*.vcf.gz'), path('*.vcf.gz.tbi')
    script:
    """
    wget $vcf_link
    wget $idx_link
    """

}

process subsetVCF {
    clusterOptions '--mem=2G --nodes=1 --cpus-per-task=1 --ntasks=1 --time=00:05:00'
    publishDir "${params.output_dir_vi}", mode: 'symlink'
    input:
        tuple path(full_vcf), path(full_idx)
    output:
        tuple path('*.vcf.gz'), path('*.vcf.gz.tbi')

    script:
    // chr_string = params.chr.collect{"chr" + it}.join("|")
    chr_string = "chr${params.chromosomes}"
    // chr_name_match = params.chr.collect{"chr" + it + " " + it}.join("\n")
    chr_name_match = "chr${params.chromosomes} ${params.chromosomes}"
    """
    echo "${chr_name_match}" > chr_name_match.txt
    module load bcftools
    # bcftools view -r "$chr_string" -Oz $full_vcf | bcftools annotate -Oz --rename-chrs chr_name_match.txt -o "${params.chromosomes.join(".")}_${full_vcf}"
    bcftools view -r "$chr_string" -Oz $full_vcf | bcftools annotate -Oz --rename-chrs chr_name_match.txt -o "${params.chromosomes}_${full_vcf}"
    # bcftools index -t "${params.chromosomes.join(".")}_${full_vcf}"
    bcftools index -t "${params.chromosomes}_${full_vcf}"
    """
}


process getConsensus {
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

workflow VARIANT_INPUT{
    take:
        input_ref
        input_vcf
        input_vcf_idx
    main:
        // Download reference file(s) for each chromosome and benchmark vcfs + index 

        downloaded_ref = getRef(input_ref, params.chromosomes)

        // Combine reference files into one file if necessary
        // if (params.chromosomes.size() > 1){
        //     ref = combineRef(input_ref.collect())
        // }
        // else{
        //     ref = downloaded_ref
        // }
        ref = downloaded_ref

        // Download benchmark vcf set and index
        downloaded_vcf_idx = getVCF(input_vcf, input_vcf_idx)


        // Subset VCF to just the chromosomes that we are sequencing and index subsetted vcf
        subset_vcf = subsetVCF(downloaded_vcf_idx)

        // Create a consensus sequence for each 'individual'
        ind_genome = getConsensus(subset_vcf, downloaded_ref.ref_path)
        
    emit:
        ref
        subsetVCF
        ind_genome
} 

workflow {
    input_ref = Channel.of(params.link_ref)
    input_vcf = Channel.of(params.link_vcf)
    input_vcf_idx = Channel.of(params.link_vcf_idx)

    VARIANT_INPUT(input_ref, input_vcf, input_vcf_idx)
}