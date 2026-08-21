#!/usr/bin/env nextflow

include {HAPPY_HAPPY} from "../../modules/nf-core/happy/happy/main.nf"
include {HAPPY_PREPY} from "../../modules/nf-core/happy/prepy/main.nf"
include {CLAIR3} from "../../modules/nf-core/clair3/main.nf"
include {SNIFFLES} from "../../modules/nf-core/sniffles/main.nf"
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main.nf'
include { getFile } from '../../modules/local/getFile.nf'
include { getFile as getFile_strat } from '../../modules/local/getFile.nf'
include { extract_bc_meta } from '../../modules/local/extractbcmeta.nf'
include { extract_bc_meta as extract_bc_meta_vcf} from '../../modules/local/extractbcmeta.nf'
include { extractTar } from '../../modules/local/extract_tar_get_strat.nf'

workflow BENCHMARK_HUMAN_VCF{
    take:
        merged_bam    // tuple (meta, bam, bai)
        ref             // path reference
        truth_vcf
    main:

        // Model and platform presets   
        model_presets = [
            ont: 'ont',
            pacbio: 'hifi',
            hifi: 'hifi'
        ]

        platform_presets = [
            ont: 'ont',
            pacbio: 'hifi',
            hifi: 'hifi'
        ]

        // Resolve model and platform
        resolved_model = params.clair3_model ?: model_presets[params.sequencing_platform] ?: 'ont'
        resolved_platform = platform_presets[params.sequencing_platform] ?: 'ont'

        // Prepare input channel with resolved values
        ch_input_clair3 = merged_bam
            .map { meta, bam, bai ->
                tuple(
                    meta,
                    bam,
                    bai,
                    resolved_model,      // packaged_model
                    [],                  // user_model (empty for packaged models)
                    resolved_platform    // platform
                )
            }


        fasta = ref.map{ it -> tuple([id:it.simpleName], it)}.collect() 
        fai_in_ch = fasta.map{meta, fa -> tuple(meta, fa, [])}
        fai = SAMTOOLS_FAIDX(fai_in_ch, false)
        fai_ch = fai.fai.collect()

        snp_calls = CLAIR3(
            ch_input_clair3,    // tuple(meta, bam, bai, packaged_model, user_model, platform)
            fasta,              // tuple(meta2, fasta)
            fai_ch          // tuple(meta3, fai)
        )

        // sv_calls = SNIFFLES
        bed_ch = channel.of(params.link_bed)
        bed_files = getFile(bed_ch.flatten())
        bed = bed_files.map{
            bed ->         
            def labels = params.barcodes.keySet().join("|")
            def matcher = bed =~ "$labels"
            def long_bc = params.barcodes[matcher[0]]
            def subs = long_bc.substring(7) as Integer as String
            def meta = [bc: "bc".plus(subs)]
            tuple(meta, bed)}

        happrep_in_vcf = snp_calls.vcf
        .map{meta, vcf -> tuple([bc:vcf.simpleName.split(/_/)[2]], meta, vcf)}
        .combine(bed, by:0)
        .map{meta1, meta2, vcf, bd -> tuple(meta2+meta1,vcf,bd)}

        prepped = HAPPY_PREPY(happrep_in_vcf, fasta, fai_ch)

        // Prep input for HAPPY_HAPPY (variant benchmark)
        truth_vcf_bc = truth_vcf.map{
            vcf ->         
            def labels = params.barcodes.keySet().join("|")
            def matcher = vcf =~ "$labels"
            def long_bc = params.barcodes[matcher[0]]
            def subs = long_bc.substring(7) as Integer as String
            def meta = [bc: "bc".plus(subs)]
            tuple(meta, vcf)}

        happy_input = prepped.preprocessed_vcf
            .map{meta, vcf -> tuple([bc:vcf.simpleName.split(/_/)[2]], meta, vcf)}
            .combine(truth_vcf_bc, by:0)
            .combine(bed, by:0)
            .map{_bc, meta, query_vcf, true_vcf, bd -> tuple(meta,query_vcf, true_vcf[0],[], bd)}

        // // No prepy
        // happy_input = snp_calls.vcf
        //     .map{meta, vcf -> tuple([bc:vcf.simpleName.split(/_/)[2]], meta, vcf)}
        //     .combine(truth_vcf_bc, by:0)
        //     .combine(bed, by:0)
        //     .map{_bc, meta, query_vcf, true_vcf, bd -> tuple(meta,query_vcf, true_vcf[0],[], bd)}

        if (file(params.stratification_tar).exists()){
            strat_dir_tar = channel.of(params.stratification_tar)
        } else {
            strat_link_ch = channel.of(params.stratification_tar)
            strat_dir_tar = getFile_strat(strat_link_ch)
        }
        
        strat_file_ch = channel.of(params.stratification_file)
        strat = extractTar(strat_dir_tar, strat_file_ch)
        strat_dir = strat.path.map{dir -> tuple([id:dir.simpleName], dir)}
        strat_tsv = strat.strat.map{tsv -> tuple([id:tsv.simpleName], tsv)}


        result = HAPPY_HAPPY(happy_input, fasta.collect(), fai_ch.collect(), [[:], []].collect() ,strat_tsv.collect(), strat_dir.collect())
        res = result.summary_csv.map{_meta, csv -> csv}
    emit:
        res
    
}