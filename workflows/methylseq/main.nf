/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap              } from 'plugin/nf-schema'
include { FASTQC                        } from '../../modules/nf-core/fastqc/main'
include { TRIMGALORE                    } from '../../modules/nf-core/trimgalore/main'
include { QUALIMAP_BAMQC                } from '../../modules/nf-core/qualimap/bamqc/main'
include { PRESEQ_LCEXTRAP               } from '../../modules/nf-core/preseq/lcextrap/main'
include { MULTIQC                       } from '../../modules/nf-core/multiqc/main'
include { CAT_FASTQ                     } from '../../modules/nf-core/cat/fastq/main'
include { FASTQ_ALIGN_DEDUP_BISMARK     } from '../../subworkflows/nf-core/fastq_align_dedup_bismark/main'
include { FASTQ_ALIGN_DEDUP_BWAMETH     } from '../../subworkflows/nf-core/fastq_align_dedup_bwameth/main'
include { FASTQ_ALIGN_DEDUP_BWAMEM      } from '../../subworkflows/nf-core/fastq_align_dedup_bwamem/main'
include { PICARD_MARKDUPLICATES         } from '../../modules/nf-core/picard/markduplicates/main'
include { PICARD_ADDORREPLACEREADGROUPS } from '../../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/picard/createsequencedictionary/main'
include { SAMTOOLS_INDEX                } from '../../modules/nf-core/samtools/index/main'
include { paramsSummaryMultiqc          } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML        } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText        } from '../../subworkflows/local/utils_nfcore_methylseq_pipeline'
include { validateInputSamplesheet      } from '../../subworkflows/local/utils_nfcore_methylseq_pipeline'
include { BAM_TAPS_CONVERSION           } from '../../subworkflows/nf-core/bam_taps_conversion'
include { BAM_METHYLDACKEL              } from '../../subworkflows/nf-core/bam_methyldackel/main'
include { TARGETED_SEQUENCING           } from '../../subworkflows/local/targeted_sequencing'

// added for SNP calling
include { SAMTOOLS_CALMD                } from '../../modules/nf-core/samtools/calmd/main'
include { REVELIO                       } from '../../modules/local/revelio/main'
include { SAMTOOLS_FAIDX                }  from '../../modules/nf-core/samtools/faidx/main'
include { GATK4_HAPLOTYPECALLER          }  from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOMICSDBIMPORT         }  from '../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS           }  from '../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATK4_MERGEVCFS               } from '../../modules/nf-core/gatk4/mergevcfs/main'
include { BCFTOOLS_SORT                 } from '../../modules/nf-core/bcftools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_GATK } from '../../modules/nf-core/samtools/index/main'
include { FQ_SUBSAMPLE                  } from '../../modules/nf-core/fq/subsample/main'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METHYLSEQ {

    take:
    samplesheet        // channel: [ path(samplesheet.csv) ]
    ch_versions        // channel: [ path(versions.yml)    ]
    ch_fasta           // channel: [ path(fasta)           ]
    ch_fasta_index     // channel: [ path(fasta index)     ]
    ch_bismark_index   // channel: [ path(bismark index)   ]
    ch_bwameth_index   // channel: [ path(bwameth index)   ]
    ch_bwamem_index    // channel: [ path(bwamem_index)    ]

    main:
    ch_fastq         = channel.empty()
    ch_fastqc_html   = channel.empty()
    ch_fastqc_zip    = channel.empty()
    ch_reads         = channel.empty()
    ch_bam           = channel.empty()
    ch_bai           = channel.empty()
    ch_gzi           = channel.empty()
    ch_bedgraph      = channel.empty()
    ch_aligner_mqc   = channel.empty()
    ch_rastair_mbias = channel.empty()
    ch_rastair_call  = channel.empty()
    ch_methylkit     = channel.empty()
    ch_mbias         = channel.empty()
    ch_qualimap      = channel.empty()
    ch_preseq        = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_bam_calmd     = channel.empty()
    ch_bai_calmd     = channel.empty()
    ch_bam_masked    = channel.empty()
    ch_bai_masked    = channel.empty()
    ch_fasta_dict    = channel.empty()
    ch_sample_vcf    = channel.empty()
    ch_gendb_input   = channel.empty()

    def callers     = params.variant_calling ? params.variant_calling.split(',')*.trim() as Set : [] as Set
    def run_masking = 'gatk' in callers || params.save_masked_bam
    // def run_masking = params.variant_calling || params.save_masked_bam
    def run_subsampling = params.subsample_fastqs != null && params.subsample_fastqs > 0
    def excluded = params.exclude_contigs ? params.exclude_contigs.split(',') as List : []


    //
    // Branch channels from input samplesheet channel
    //
    ch_samplesheet = samplesheet
                        .branch { meta, fastqs ->
                            single  : fastqs.size() == 1
                                return [ meta, fastqs.flatten() ]
                            multiple: fastqs.size() > 1
                                return [ meta, fastqs.flatten() ]
                        }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_samplesheet.multiple
    )
    ch_fastq    = CAT_FASTQ.out.reads.mix(ch_samplesheet.single)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    if (run_subsampling) {

        FQ_SUBSAMPLE ( ch_fastq )

        ch_fastq = FQ_SUBSAMPLE.out.fastq

    }

    //
    // MODULE: Run FastQC
    //
    if (!params.skip_fastqc) {
        FASTQC (
            ch_fastq
        )
        ch_fastqc_html   = FASTQC.out.html
        ch_fastqc_zip    = FASTQC.out.zip
    } else {
        ch_fastqc_html   = channel.empty()
        ch_fastqc_zip    = channel.empty()
    }

    //
    // MODULE: Run TrimGalore!
    //
    if (!params.skip_trimming) {
        TRIMGALORE(
            ch_fastq
        )
        ch_reads    = TRIMGALORE.out.reads
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
    } else {
        ch_reads    = ch_fastq
    }

    //
    // SUBWORKFLOW: Align reads, deduplicate and extract methylation with Bismark
    //

    if (params.taps && params.aligner != 'bwamem') {
        log.info "TAPS protocol detected and aligner is not 'bwamem'. We recommend using bwa-mem for TAPS protocol as it is optimized for this type of data."
    }

    // Aligner: bismark or bismark_hisat
    if (params.aligner =~ /bismark/ ) {
        //
        // Run Bismark alignment + downstream processing
        //
        ch_bismark_inputs = ch_reads
            .combine(ch_fasta)
            .combine(ch_bismark_index)
            .multiMap { meta, reads, meta_fasta, fasta, meta_bismark, bismark_index ->
                reads: [ meta, reads ]
                fasta: [ meta_fasta, fasta ]
                bismark_index: [ meta_bismark, bismark_index ]
            }

        FASTQ_ALIGN_DEDUP_BISMARK (
            ch_bismark_inputs.reads,
            ch_bismark_inputs.fasta,
            ch_bismark_inputs.bismark_index,
            params.skip_deduplication || params.rrbs,
            params.cytosine_report || params.nomeseq
        )
        ch_bam         = FASTQ_ALIGN_DEDUP_BISMARK.out.bam
        ch_bai         = FASTQ_ALIGN_DEDUP_BISMARK.out.bai
        ch_bedgraph    = FASTQ_ALIGN_DEDUP_BISMARK.out.methylation_bedgraph
        ch_aligner_mqc = FASTQ_ALIGN_DEDUP_BISMARK.out.multiqc
        ch_versions    = ch_versions.mix(FASTQ_ALIGN_DEDUP_BISMARK.out.versions)
    }
    // Aligner: bwameth
    else if (params.aligner == 'bwameth'){

        ch_bwameth_inputs = ch_reads
            .combine(ch_fasta)
            .combine(ch_fasta_index)
            .combine(ch_bwameth_index)
            .multiMap { meta, reads, meta_fasta, fasta, meta_fasta_index, fasta_index, meta_bwameth, bwameth_index ->
                reads: [ meta, reads ]
                fasta: [ meta_fasta, fasta ]
                fasta_index: [ meta_fasta_index, fasta_index ]
                bwameth_index: [ meta_bwameth, bwameth_index ]
            }

        FASTQ_ALIGN_DEDUP_BWAMETH (
            ch_bwameth_inputs.reads,
            ch_bwameth_inputs.fasta,
            ch_bwameth_inputs.fasta_index,
            ch_bwameth_inputs.bwameth_index,
            params.skip_deduplication || params.rrbs,
            workflow.profile.tokenize(',').intersect(['gpu']).size() >= 1
        )
        ch_bam         = FASTQ_ALIGN_DEDUP_BWAMETH.out.bam
        ch_bai         = FASTQ_ALIGN_DEDUP_BWAMETH.out.bai
        ch_aligner_mqc = FASTQ_ALIGN_DEDUP_BWAMETH.out.multiqc
        ch_versions    = ch_versions.mix(FASTQ_ALIGN_DEDUP_BWAMETH.out.versions)
    }

    // Aligner: bwamem
    else if (params.aligner == 'bwamem'){

        ch_bwamem_inputs = ch_reads
            .combine(ch_fasta)
            .combine(ch_fasta_index)
            .combine(ch_bwamem_index)
            .multiMap { meta, reads, meta_fasta, fasta, meta_fasta_index, fasta_index, meta_bwamem, bwamem_index ->
                reads: [ meta, reads ]
                fasta: [ meta_fasta, fasta ]
                fasta_index: [ meta_fasta_index, fasta_index ]
                bwamem_index: [ meta_bwamem, bwamem_index ]
            }

        FASTQ_ALIGN_DEDUP_BWAMEM (
            ch_bwamem_inputs.reads,
            ch_bwamem_inputs.fasta,
            ch_bwamem_inputs.fasta_index,
            ch_bwamem_inputs.bwamem_index,
            params.skip_deduplication
        )

        ch_bam         = FASTQ_ALIGN_DEDUP_BWAMEM.out.bam
        ch_bai         = FASTQ_ALIGN_DEDUP_BWAMEM.out.bai
        ch_aligner_mqc = FASTQ_ALIGN_DEDUP_BWAMEM.out.multiqc
        ch_versions    = ch_versions.mix(FASTQ_ALIGN_DEDUP_BWAMEM.out.versions.unique{ it.baseName })
    }

    else {
        error "ERROR: Invalid aligner '${params.aligner}'. Valid options are: 'bismark', 'bismark_hisat', 'bwameth' or 'bwamem'."
    }

    //
    // Subworkflow: Count positive mC->T conversion rates as a readout for DNA methylation
    //
    if (params.taps || params.aligner == 'bwamem') {

        ch_bam_bai = ch_bam.join(ch_bai)
        ch_taps_inputs = ch_bam_bai
            .combine(ch_fasta)         // broadcast fasta
            .combine(ch_fasta_index)   // broadcast fai
            .multiMap { meta, bam, bai, _meta_fasta, fasta, _meta_fai, fai ->
                bam:         [ meta, bam ]      // use sample meta so subworkflow aligns properly
                bai:         [ meta, bai ]      // use sample meta so subworkflow aligns properly
                fasta:       [ meta, fasta ]    // use sample meta so subworkflow aligns properly
                fasta_index: [ meta, fai ]      // same here
            }

        BAM_TAPS_CONVERSION(
            ch_taps_inputs.bam,
            ch_taps_inputs.bai,
            ch_taps_inputs.fasta,
            ch_taps_inputs.fasta_index
        )
        ch_rastair_mbias = BAM_TAPS_CONVERSION.out.mbias // channel: [ val(meta), [ txt ] ]
        ch_rastair_call  = BAM_TAPS_CONVERSION.out.call // channel: [ val(meta), [ txt ] ]
        ch_versions      = ch_versions.mix(BAM_TAPS_CONVERSION.out.versions)
    }

    //
    // Subworkflow: Count negative C->T conversion rates as a readout for DNA methylation
    //
    else if (!params.taps && params.aligner == 'bwameth') {

        ch_bam_bai = ch_bam.join(ch_bai)
        ch_methyldackel_inputs = ch_bam_bai
            .combine(ch_fasta)         // broadcast fasta
            .combine(ch_fasta_index)   // broadcast fai
            .multiMap { meta, bam, bai, _meta_fasta, fasta, _meta_fai, fai ->
                bam:         [ meta, bam ]      // use sample meta so subworkflow aligns properly
                bai:         [ meta, bai ]      // use sample meta so subworkflow aligns properly
                fasta:       [ meta, fasta ]    // use sample meta so subworkflow aligns properly
                fasta_index: [ meta, fai ]      // same here
            }

        BAM_METHYLDACKEL (
            ch_methyldackel_inputs.bam,
            ch_methyldackel_inputs.bai,
            ch_methyldackel_inputs.fasta,
            ch_methyldackel_inputs.fasta_index
        )
        ch_bedgraph    = BAM_METHYLDACKEL.out.methydackel_extract_bedgraph  // channel: [ val(meta), [ bedgraph ]  ]
        ch_methylkit   = BAM_METHYLDACKEL.out.methydackel_extract_methylkit // channel: [ val(meta), [ methylkit ] ]
        ch_mbias       = BAM_METHYLDACKEL.out.methydackel_mbias // channel: [ val(meta), [ mbias ] ]
        ch_versions    = ch_versions.mix(BAM_METHYLDACKEL.out.versions)
    }

    //
    // MODULE: Qualimap BamQC
    // skipped by default. to use run with `--run_qualimap` param.
    //
    if(params.run_qualimap) {
        QUALIMAP_BAMQC (
            ch_bam,
            params.bamqc_regions_file ? channel.fromPath( params.bamqc_regions_file, checkIfExists: true ).toList() : []
        )
        ch_qualimap = QUALIMAP_BAMQC.out.results
        ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions)
    }

    //
    // MODULE: Targeted sequencing analysis
    // skipped by default. to use run with `--run_targeted_sequencing` param.
    //
    if (params.run_targeted_sequencing){
        if (params.taps || params.aligner == 'bwamem') {
            error "ERROR: --run_targeted_sequencing can't be running using rastair (methylation caller for TAPS) "
        }
        if (!params.target_regions_file) {
            error "ERROR: --target_regions_file must be specified when using --run_targeted_sequencing"
        }
        TARGETED_SEQUENCING (
            ch_bedgraph,
            channel.fromPath(params.target_regions_file, checkIfExists: true),
            ch_fasta,
            ch_fasta_index,
            ch_bam,
            ch_bai,
            ch_gzi,
            params.collecthsmetrics
        )
        ch_versions = ch_versions.mix(TARGETED_SEQUENCING.out.versions)
    }

    //
    // MODULE: Preseq LCEXTRAP
    // skipped by default. to use run with `--run_preseq` param.
    //
    if(params.run_preseq) {
        PRESEQ_LCEXTRAP (
            ch_bam
        )
        ch_preseq   = PRESEQ_LCEXTRAP.out.lc_extrap
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'methylseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // Topic channel versions - written separately to avoid blocking MULTIQC
    // These will be merged into the main versions file on workflow completion
    //
    channel.topic("versions")
        .distinct()
        .filter { entry -> !(entry instanceof Path) }
        .map { process, tool, version ->
            def processName = process[process.lastIndexOf(':')+1..-1]
            "${processName}:\n  ${tool}: ${version}"
        }
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_methylseq_topic_versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: samtools calmd
    //
    if (run_masking) {

        SAMTOOLS_FAIDX (
            ch_fasta,
            [[:], []],
            false
        )

        ch_fasta_index = SAMTOOLS_FAIDX.out.fai

        ch_calmd_inputs = ch_bam
            .combine(ch_fasta)
            .combine(ch_fasta_index)
            .multiMap { meta, bam, _meta_fasta, fasta, _meta_fai, fai ->
                bam:   [ meta, bam ]
                fasta: [ meta, fasta, fai ]
            }
        
        SAMTOOLS_CALMD (
            ch_calmd_inputs.bam,
            ch_calmd_inputs.fasta
        )

        ch_bam_calmd = SAMTOOLS_CALMD.out.bam

    }


    //
    // MODULE: samtools index
    //
    if (run_masking) {
        
        SAMTOOLS_INDEX (
            ch_bam_calmd
        )

        ch_bai_calmd = SAMTOOLS_INDEX.out.bai

    }

    //
    // MODULE: revelio
    //
    if (run_masking) {
        
        ch_revelio_input = ch_bam_calmd.join(ch_bai_calmd)

        REVELIO (
            ch_revelio_input
        )

        ch_bam_masked = REVELIO.out.masked_bam

    }


    //
    // SNP Calling with GATK
    //
    if ('gatk' in callers) {
        

        // requires read groups adding to revelio output
        PICARD_ADDORREPLACEREADGROUPS (
            ch_bam_masked,
            ch_fasta.first(),       // fasta (required by module signature)
            ch_fasta_index.first()  // fai
        )

        // index masked bams with RGs
        SAMTOOLS_INDEX_GATK (
            PICARD_ADDORREPLACEREADGROUPS.out.bam
        )

        ch_bai_masked = SAMTOOLS_INDEX_GATK.out.bai
        
        PICARD_CREATESEQUENCEDICTIONARY (
            ch_fasta
        )

        ch_fasta_dict = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict

        ch_gatk_inputs = PICARD_ADDORREPLACEREADGROUPS.out.bam
            .join(ch_bai_masked)
            .combine(ch_fasta)
            .combine(ch_fasta_index)
            .combine(ch_fasta_dict)
            .multiMap { meta, bam_masked, bai_masked, _meta_fasta, fasta, _meta_fai, fai, _meta_dict, dict ->
                bam:   [ meta, bam_masked, bai_masked, [], [] ]
                fasta: [ _meta_fasta, fasta ]
                fai:   [ _meta_fai, fai ]
                dict:  [ _meta_dict, dict ]
            }
        
        //
        // MODULE: gatk haplotypecaller
        //
        
        GATK4_HAPLOTYPECALLER (
            ch_gatk_inputs.bam,   // [ meta, bam, bai, [], [] ]
            ch_gatk_inputs.fasta, // [ meta2, fasta ]
            ch_gatk_inputs.fai,   // [ meta3, fai ]
            ch_gatk_inputs.dict,  // [ meta4, dict ]
            [ [:], [] ],          // dbsnp
            [ [:], [] ]           // dbsnp_tbi
        )

        ch_sample_vcf = GATK4_HAPLOTYPECALLER.out.vcf


        //
        // MODULE: gatk genomics DB import
        //
        ch_intervals = ch_fasta_index
            .map { meta, fai -> fai }
            .splitCsv(sep: '\t')
            .map { row -> row[0] }
            .filter { !(it in excluded) }
            .collectFile(name: 'genome.intervals', newLine: true)

        ch_gendb_input = GATK4_HAPLOTYPECALLER.out.vcf
            .join(GATK4_HAPLOTYPECALLER.out.tbi)
            .map { meta, vcf, tbi -> [ vcf, tbi ] }
            .collect()
            .map { files ->
                def vcfs = files.collate(2).collect { it[0] }
                def tbis = files.collate(2).collect { it[1] }
                [ [id: 'joint'], vcfs, tbis ]
            }
            .combine(ch_intervals)
            .map { meta, vcfs, tbis, intervals ->
                [ meta, vcfs, tbis, intervals, [], [] ]
            }

        GATK4_GENOMICSDBIMPORT(
                ch_gendb_input,
                false,   // run_intlist
                false,   // run_updatewspace
                false    // input_map
            )

        genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{ 
            meta, genomicsdb -> [ meta, genomicsdb, [], [], [] ] // add three empty slots so output is ready for genotype gVCFs
             }

        //
        // MODULE: gatk genotype gvcfs
        //

        // Joint genotyping performed using GenotypeGVCFs

           GATK4_GENOTYPEGVCFS(
                genotype_input,       // [ meta, genomicsdb, [], [], [] ]
                ch_fasta,             // [ meta2, fasta ]
                ch_fasta_index,       // [ meta3, fai ]
                ch_fasta_dict,        // [ meta4, dict ]
                [ [:], [] ],          // dbsnp
                [ [:], [] ]           // dbsnp_tbi
            )


        // Sort and merge joint genotype calls 

        // BCFTOOLS_SORT(GATK4_GENOTYPEGVCFS.out.vcf)
        // ch_joint_vcf  = BCFTOOLS_SORT.out.vcf

        // Merge scatter/gather vcfs & index
        // Rework meta for variantscalled.csv and annotation tools
        // GATK4_MERGEVCFS(GATK4_GENOTYPEGVCFS.out.vcf, ch_fasta_dict)





    } // end if gatk in callers



    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        ch_multiqc_config        = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ?
            channel.fromPath(params.multiqc_config, checkIfExists: true) :
            channel.empty()
        ch_multiqc_logo          = params.multiqc_logo ?
            channel.fromPath(params.multiqc_logo, checkIfExists: true) :
            channel.empty()

        summary_params           = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary      = channel.value(paramsSummaryMultiqc(summary_params))

        ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
            file(params.multiqc_methods_description, checkIfExists: true) :
            file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_methods_description.collectFile(
                name: 'methods_description_mqc.yaml',
                sort: true
            )
        )

        if(params.run_qualimap) {
            ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.collect{ it[1] }.ifEmpty([]))
        }
        if (params.run_preseq) {
            ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.log.collect{ it[1] }.ifEmpty([]))
        }
        ch_multiqc_files = ch_multiqc_files.mix(ch_aligner_mqc.ifEmpty([]))
        if (!params.skip_trimming) {
            ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{ it[1] })
        }
        if (params.run_targeted_sequencing) {
            if (params.collecthsmetrics) {
                ch_multiqc_files = ch_multiqc_files.mix(TARGETED_SEQUENCING.out.picard_metrics.collect{ it[1] }.ifEmpty([]))
            }
        }
        if (!params.skip_fastqc) {
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ it[1] }.ifEmpty([]))
        }

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )
        ch_multiqc_report = MULTIQC.out.report.toList()
    } else {
        ch_multiqc_report = channel.empty()
    }


    emit:
    bam            = ch_bam                      // channel: [ val(meta), path(bam) ]
    bai            = ch_bai                      // channel: [ val(meta), path(bai) ]
    rastair_mbias  = ch_rastair_mbias            // channel: [ val(meta), path(rastair_mbias) ]
    rastair_call   = ch_rastair_call             // channel: [ val(meta), path(rastair_call) ]
    methylkit      = ch_methylkit                 // channel: [ val(meta), path(methylkit) ]
    mbias          = ch_mbias                     // channel: [ val(meta), path(mbias) ]
    bedgraph       = ch_bedgraph                  // channel: [ val(meta), path(bedgraph) ]
    qualimap       = ch_qualimap                 // channel: [ val(meta), path(qualimap) ]
    preseq         = ch_preseq                   // channel: [ val(meta), path(preseq) ]
    multiqc_report = ch_multiqc_report            // channel: [ path(multiqc_report.html )  ]
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    bam_calmd      = ch_bam_calmd                // channel: [ val(meta), path(bam) ]
    bai_calmd      = ch_bai_calmd                // channel: [ val(meta), path(bai) ]
    bam_masked     = ch_bam_masked                // channel: [ val(meta), path(bam) ]
    bai_masked     = ch_bai_masked               // channel: [ val(meta), path(bai) ]
    fasta_dict     = ch_fasta_dict               // channel: [ val(meta2), path(dict) ]
    sample_vcf     = ch_sample_vcf               // channel: [ val(meta), path(vcf.gz) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
