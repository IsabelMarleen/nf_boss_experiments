#!/usr/bin/env nextflow

// Run using 'nextflow run boss_variant_input.nf -entry VARIANT_INPUT -with-conda -resume -c ../nextflow.config'

/*
 * Pipeline parameters in nextflow.config
 */


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
    ref_base = "${refs.first().getBaseName(3)}"
    """
    cat ${refs} > "${ref_base}.${params.chromosomes.join(".")}.fa.gz"
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
    if (params.chromosomes.size() > 1){
        chr_string = params.chromosomes.collect{"chr" + it}.join(",")
        chr_name_match = params.chromosomes.collect{"chr" + it + " " + it}.join("\n")
        chr_prefix = params.chromosomes.join("_")
    }
    else{
        chr_string = "chr${params.chromosomes.first()}"
        chr_name_match = "chr${params.chromosomes.first()} ${params.chromosomes.first()}"
        chr_prefix = params.chromosomes.first()
    }
    """
    echo "${chr_name_match}" > chr_name_match.txt
    module load bcftools
    bcftools view -r "$chr_string" -Oz $full_vcf | bcftools annotate -Oz --rename-chrs chr_name_match.txt -o "${chr_prefix}_${full_vcf}"
    bcftools index -t "${chr_prefix}_${full_vcf}"
    """
}


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

workflow VARIANT_INPUT{
    take:
        input_ref
        input_vcf
        input_vcf_idx
    main:
        // Download reference file(s) for each chromosome and benchmark vcfs + index 
        chrs = Channel.from(params.chromosomes)
        downloaded_ref = getRef(input_ref.first(), chrs)

        // Combine reference files into one file if necessary
        if (params.chromosomes.size() > 1){
            ref = combineRef(downloaded_ref.collect())
        }
        else{
            ref = downloaded_ref
        }

        // Download benchmark vcf set and index
        downloaded_vcf_idx = getVCF(input_vcf, input_vcf_idx)


        // Subset VCF to just the chromosomes that we are sequencing and index subsetted vcf
        subset_vcf = subsetVCF(downloaded_vcf_idx)

        // Create a consensus sequence for each 'individual'
        ind_genome = getConsensus(subset_vcf, ref)
        
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