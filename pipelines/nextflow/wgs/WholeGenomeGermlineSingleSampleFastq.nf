#!/usr/bin/env nextflow

/*
 * Nextflow implementation of the WARP WholeGenomeGermlineSingleSampleFastq pipeline
 * 
 * This pipeline implements data pre-processing and initial variant calling (GVCF
 * generation) according to the GATK Best Practices for germline SNP and
 * Indel discovery in human whole-genome data.
 *
 * Original WDL: https://github.com/broadinstitute/warp/blob/master/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSampleFastq.wdl
 *
 * Requirements/expectations:
 * - Human whole-genome pair-end sequencing data in FASTQ format
 * - GVCF output names must end in ".g.vcf.gz"
 * - Reference genome must be Hg38 with ALT contigs
 */

nextflow.enable.dsl = 2

// Pipeline version
def pipeline_version = "1.0.0-nf"

/*
 * Pipeline parameters with defaults
 */
params {
    // Input files
    input_fastq_r1       = null
    input_fastq_r2       = null
    sample_name          = null
    base_file_name       = null
    final_gvcf_base_name = null
    
    // Read group information
    read_group_id        = null
    read_group_platform  = "ILLUMINA"
    read_group_pu        = null
    read_group_library   = null
    read_group_center    = null
    
    // Reference files
    reference_fasta      = null
    reference_fasta_index = null
    reference_dict       = null
    reference_alt        = null
    reference_amb        = null
    reference_ann        = null
    reference_bwt        = null
    reference_pac        = null
    reference_sa         = null
    
    // Known sites
    dbsnp_vcf            = null
    dbsnp_vcf_index      = null
    known_indels_sites_vcfs = []
    known_indels_sites_indices = []
    
    // Interval lists
    calling_interval_list = null
    evaluation_interval_list = null
    wgs_coverage_interval_list = null
    
    // Contamination files
    contamination_sites_ud = null
    contamination_sites_bed = null
    contamination_sites_mu = null
    haplotype_database_file = null
    
    // Optional files
    fingerprint_genotypes_file = null
    fingerprint_genotypes_index = null
    
    // Optional DRAGMAP reference
    dragmap_reference_bin = null
    dragmap_hash_table_cfg_bin = null
    dragmap_hash_table_cmp = null
    
    // Pipeline options
    provide_bam_output = false
    use_gatk3_haplotype_caller = true
    dragen_functional_equivalence_mode = false
    dragen_maximum_quality_mode = false
    run_dragen_mode_variant_calling = false
    use_spanning_event_genotyping = true
    unmap_contaminant_reads = true
    perform_bqsr = true
    use_bwa_mem = true
    use_dragen_hard_filtering = false
    
    // Scatter settings
    haplotype_scatter_count = 50
    break_bands_at_multiples_of = 1000000
    
    // Output directory
    outdir = "results"
    
    // Resource settings
    max_cpus = 16
    max_memory = "64.GB"
    max_time = "240.h"
}

/*
 * Validate required parameters
 */
if (!params.input_fastq_r1) error "Missing required parameter: input_fastq_r1"
if (!params.input_fastq_r2) error "Missing required parameter: input_fastq_r2"
if (!params.sample_name) error "Missing required parameter: sample_name"
if (!params.reference_fasta) error "Missing required parameter: reference_fasta"

/*
 * Set derived parameters
 */
def base_file_name = params.base_file_name ?: params.sample_name
def final_gvcf_base_name = params.final_gvcf_base_name ?: base_file_name
def read_length = 250
def lod_threshold = -20.0
def cross_check_fingerprints_by = "READGROUP"
def recalibrated_bam_basename = "${base_file_name}.aligned.duplicates_marked.recalibrated"

// Set DRAGEN-related arguments according to the preset arguments
def run_dragen_mode_variant_calling_ = (params.dragen_functional_equivalence_mode || params.dragen_maximum_quality_mode) ? true : params.run_dragen_mode_variant_calling
def use_spanning_event_genotyping_ = params.dragen_functional_equivalence_mode ? false : (params.dragen_maximum_quality_mode ? true : params.use_spanning_event_genotyping)
def unmap_contaminant_reads_ = params.dragen_functional_equivalence_mode ? false : (params.dragen_maximum_quality_mode ? true : params.unmap_contaminant_reads)
def perform_bqsr_ = (params.dragen_functional_equivalence_mode || params.dragen_maximum_quality_mode) ? false : params.perform_bqsr
def use_bwa_mem_ = (params.dragen_functional_equivalence_mode || params.dragen_maximum_quality_mode) ? false : params.use_bwa_mem
def use_gatk3_haplotype_caller_ = (params.dragen_functional_equivalence_mode || params.dragen_maximum_quality_mode) ? false : params.use_gatk3_haplotype_caller

/*
 * Include processes
 */
include { BWA_MEM_ALIGN } from './modules/alignment.nf'
include { MARK_DUPLICATES } from './modules/processing.nf'
include { SORT_BAM } from './modules/processing.nf'
include { BASE_RECALIBRATOR } from './modules/gatk.nf'
include { APPLY_BQSR } from './modules/gatk.nf'
include { GATHER_BAM_FILES } from './modules/processing.nf'
include { COLLECT_WGS_METRICS } from './modules/qc.nf'
include { COLLECT_RAW_WGS_METRICS } from './modules/qc.nf'
include { AGGREGATED_BAM_QC } from './modules/qc.nf'
include { BAM_TO_CRAM } from './modules/processing.nf'
include { HAPLOTYPE_CALLER } from './modules/gatk.nf'
include { CHECK_CONTAMINATION } from './modules/qc.nf'

/*
 * Main workflow
 */
workflow {
    // Validate conflicting parameters
    if (params.dragen_functional_equivalence_mode && params.dragen_maximum_quality_mode) {
        error "Both dragen_functional_equivalence_mode and dragen_maximum_quality_mode have been set to true, however, they are mutually exclusive."
    }
    
    if (run_dragen_mode_variant_calling_ && use_gatk3_haplotype_caller_) {
        error "DRAGEN mode variant calling has been activated, however, the HaplotypeCaller version has been set to use GATK 3."
    }
    
    // Create input channels
    fastq_r1_ch = Channel.fromPath(params.input_fastq_r1, checkIfExists: true)
    fastq_r2_ch = Channel.fromPath(params.input_fastq_r2, checkIfExists: true)
    
    // Reference files channels
    reference_fasta_ch = Channel.fromPath(params.reference_fasta, checkIfExists: true)
    reference_fasta_index_ch = Channel.fromPath(params.reference_fasta_index, checkIfExists: true)
    reference_dict_ch = Channel.fromPath(params.reference_dict, checkIfExists: true)
    
    // BWA index files
    reference_alt_ch = Channel.fromPath(params.reference_alt, checkIfExists: true)
    reference_amb_ch = Channel.fromPath(params.reference_amb, checkIfExists: true)
    reference_ann_ch = Channel.fromPath(params.reference_ann, checkIfExists: true)
    reference_bwt_ch = Channel.fromPath(params.reference_bwt, checkIfExists: true)
    reference_pac_ch = Channel.fromPath(params.reference_pac, checkIfExists: true)
    reference_sa_ch = Channel.fromPath(params.reference_sa, checkIfExists: true)
    
    // Combine FASTQ files
    fastq_pairs_ch = fastq_r1_ch.combine(fastq_r2_ch)
    
    // Alignment step
    BWA_MEM_ALIGN(
        fastq_pairs_ch,
        reference_fasta_ch,
        reference_fasta_index_ch,
        reference_alt_ch,
        reference_amb_ch,
        reference_ann_ch,
        reference_bwt_ch,
        reference_pac_ch,
        reference_sa_ch,
        params.sample_name,
        params.read_group_id ?: params.sample_name,
        params.read_group_platform,
        params.read_group_pu ?: "unknown",
        params.read_group_library ?: params.sample_name,
        params.read_group_center ?: "unknown"
    )
    
    // Mark duplicates
    MARK_DUPLICATES(
        BWA_MEM_ALIGN.out.aligned_bam,
        base_file_name
    )
    
    // Sort BAM
    SORT_BAM(
        MARK_DUPLICATES.out.marked_bam,
        base_file_name
    )
    
    // Base Quality Score Recalibration (if enabled)
    if (perform_bqsr_) {
        // Prepare known sites
        known_sites_ch = Channel.fromPath(params.dbsnp_vcf, checkIfExists: true)
            .mix(Channel.fromList(params.known_indels_sites_vcfs).flatten())
        
        known_sites_indices_ch = Channel.fromPath(params.dbsnp_vcf_index, checkIfExists: true)
            .mix(Channel.fromList(params.known_indels_sites_indices).flatten())
        
        BASE_RECALIBRATOR(
            SORT_BAM.out.sorted_bam,
            reference_fasta_ch,
            reference_fasta_index_ch,
            reference_dict_ch,
            known_sites_ch.collect(),
            known_sites_indices_ch.collect(),
            base_file_name
        )
        
        APPLY_BQSR(
            SORT_BAM.out.sorted_bam,
            BASE_RECALIBRATOR.out.recal_table,
            reference_fasta_ch,
            reference_fasta_index_ch,
            reference_dict_ch,
            base_file_name
        )
        
        final_bam_ch = APPLY_BQSR.out.recalibrated_bam
        final_bam_index_ch = APPLY_BQSR.out.recalibrated_bai
    } else {
        final_bam_ch = SORT_BAM.out.sorted_bam
        final_bam_index_ch = SORT_BAM.out.sorted_bai
    }
    
    // Quality control metrics
    wgs_interval_ch = Channel.fromPath(params.wgs_coverage_interval_list, checkIfExists: true)
    
    COLLECT_WGS_METRICS(
        final_bam_ch,
        final_bam_index_ch,
        reference_fasta_ch,
        reference_fasta_index_ch,
        wgs_interval_ch,
        base_file_name,
        read_length
    )
    
    COLLECT_RAW_WGS_METRICS(
        final_bam_ch,
        final_bam_index_ch,
        reference_fasta_ch,
        reference_fasta_index_ch,
        wgs_interval_ch,
        base_file_name,
        read_length
    )
    
    // Contamination check
    contamination_sites_ud_ch = Channel.fromPath(params.contamination_sites_ud, checkIfExists: true)
    contamination_sites_bed_ch = Channel.fromPath(params.contamination_sites_bed, checkIfExists: true)
    contamination_sites_mu_ch = Channel.fromPath(params.contamination_sites_mu, checkIfExists: true)
    
    CHECK_CONTAMINATION(
        final_bam_ch,
        final_bam_index_ch,
        contamination_sites_ud_ch,
        contamination_sites_bed_ch,
        contamination_sites_mu_ch,
        base_file_name
    )
    
    // Aggregated BAM QC
    haplotype_db_ch = Channel.fromPath(params.haplotype_database_file, checkIfExists: true)
    
    AGGREGATED_BAM_QC(
        final_bam_ch,
        final_bam_index_ch,
        reference_fasta_ch,
        reference_fasta_index_ch,
        reference_dict_ch,
        haplotype_db_ch,
        base_file_name,
        params.sample_name
    )
    
    // Convert BAM to CRAM
    BAM_TO_CRAM(
        final_bam_ch,
        final_bam_index_ch,
        reference_fasta_ch,
        reference_fasta_index_ch,
        reference_dict_ch,
        MARK_DUPLICATES.out.duplicate_metrics,
        base_file_name
    )
    
    // Variant calling
    calling_interval_ch = Channel.fromPath(params.calling_interval_list, checkIfExists: true)
    evaluation_interval_ch = Channel.fromPath(params.evaluation_interval_list, checkIfExists: true)
    dbsnp_vcf_ch = Channel.fromPath(params.dbsnp_vcf, checkIfExists: true)
    dbsnp_vcf_index_ch = Channel.fromPath(params.dbsnp_vcf_index, checkIfExists: true)
    
    HAPLOTYPE_CALLER(
        final_bam_ch,
        final_bam_index_ch,
        reference_fasta_ch,
        reference_fasta_index_ch,
        reference_dict_ch,
        calling_interval_ch,
        dbsnp_vcf_ch,
        dbsnp_vcf_index_ch,
        CHECK_CONTAMINATION.out.contamination_value,
        final_gvcf_base_name,
        use_gatk3_haplotype_caller_,
        run_dragen_mode_variant_calling_,
        use_spanning_event_genotyping_
    )
    
    // Emit outputs
    emit:
        output_bam = params.provide_bam_output ? final_bam_ch : Channel.empty()
        output_bam_index = params.provide_bam_output ? final_bam_index_ch : Channel.empty()
        output_cram = BAM_TO_CRAM.out.cram
        output_cram_index = BAM_TO_CRAM.out.cram_index
        output_cram_md5 = BAM_TO_CRAM.out.cram_md5
        output_vcf = HAPLOTYPE_CALLER.out.gvcf
        output_vcf_index = HAPLOTYPE_CALLER.out.gvcf_index
        duplicate_metrics = MARK_DUPLICATES.out.duplicate_metrics
        wgs_metrics = COLLECT_WGS_METRICS.out.wgs_metrics
        raw_wgs_metrics = COLLECT_RAW_WGS_METRICS.out.raw_wgs_metrics
        contamination_value = CHECK_CONTAMINATION.out.contamination_value
        selfSM = CHECK_CONTAMINATION.out.selfSM
        agg_alignment_summary_metrics = AGGREGATED_BAM_QC.out.alignment_summary_metrics
        agg_insert_size_metrics = AGGREGATED_BAM_QC.out.insert_size_metrics
        agg_gc_bias_metrics = AGGREGATED_BAM_QC.out.gc_bias_metrics
}

/*
 * Completion message
 */
workflow.onComplete {
    println """
    Pipeline completed!
    
    Results are located in: ${params.outdir}
    
    Main outputs:
    - GVCF: ${params.outdir}/${final_gvcf_base_name}.g.vcf.gz
    - CRAM: ${params.outdir}/${base_file_name}.cram
    - Metrics: Various QC metrics files
    """
}

/*
 * Error handling
 */
workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
