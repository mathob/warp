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

/*
 * Pipeline parameters with defaults
 */
// Pipeline parameters with defaults
params.pipeline_version = "1.0.0-nf"
// Input files
params.input_fastq_r1       = null
params.input_fastq_r2       = null
params.sample_name          = null
params.base_file_name       = null
params.final_gvcf_base_name = null

// Read group information
params.read_group_id        = null
params.read_group_platform  = "ILLUMINA"
params.read_group_pu        = null
params.read_group_library   = null
params.read_group_center    = null

// Reference files
params.reference_fasta      = null
params.reference_fasta_index = null
params.reference_dict       = null
params.reference_alt        = null
params.reference_amb        = null
params.reference_ann        = null
params.reference_bwt        = null
params.reference_pac        = null
params.reference_sa         = null

// Known sites
params.dbsnp_vcf            = null
params.dbsnp_vcf_index      = null
params.known_indels_sites_vcfs = []
params.known_indels_sites_indices = []

// Interval lists
params.calling_interval_list = null
params.evaluation_interval_list = null
params.wgs_coverage_interval_list = null

// Contamination files
params.contamination_sites_ud = null
params.contamination_sites_bed = null
params.contamination_sites_mu = null
params.haplotype_database_file = null

// Optional files
params.fingerprint_genotypes_file = null
params.fingerprint_genotypes_index = null

// Optional DRAGMAP reference
params.dragmap_reference_bin = null
params.dragmap_hash_table_cfg_bin = null
params.dragmap_hash_table_cmp = null
params.str_table_file = null

// Pipeline options
params.provide_bam_output = false
params.use_gatk3_haplotype_caller = true
params.dragen_functional_equivalence_mode = false
params.dragen_maximum_quality_mode = false
params.run_dragen_mode_variant_calling = false
params.use_spanning_event_genotyping = true
params.unmap_contaminant_reads = true
params.perform_bqsr = true
params.aligner = "dragmap"  // Options: "dragmap" (default) or "bwa-mem"
params.use_dragen_hard_filtering = false

// Scatter settings
params.haplotype_scatter_count = 50
params.break_bands_at_multiples_of = 1000000

// Output directory
params.outdir = "results"

// Resource settings
params.max_cpus = 16
params.max_memory = "64.GB"
params.max_time = "240.h"

/*
 * Include processes
 */
include { CALIBRATE_DRAGSTR_MODEL } from './modules/alignment.nf'
include { FASTQ2UBAM } from './modules/fastq2ubam.nf'
include { BWA_MEM_ALIGN } from './modules/alignment.nf'
include { DRAGMAP_ALIGN } from './modules/alignment.nf'
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
include { DRAGEN_HARD_VARIANT_FILTRATION } from './modules/gatk.nf'
include { CHECK_CONTAMINATION } from './modules/qc.nf'

/*
 * Main workflow
 */
workflow {

    /*
     * All variable definitions and logic must be inside the workflow block in DSL2.
     */

    // All variable definitions and logic must be inside the workflow block in DSL2

    def base_file_name = params.base_file_name ?: params.sample_name
    def final_gvcf_base_name = params.final_gvcf_base_name ?: base_file_name

    // Derived parameters
    def read_length = 250
    def lod_threshold = -20.0
    def cross_check_fingerprints_by = "READGROUP"
    def recalibrated_bam_basename = "${base_file_name}.aligned.duplicates_marked.recalibrated"

    // Set DRAGEN-related arguments according to the preset arguments
    def run_dragen_mode_variant_calling_ = (params.dragen_functional_equivalence_mode || params.dragen_maximum_quality_mode) ? true : params.run_dragen_mode_variant_calling
    def use_spanning_event_genotyping_ = params.dragen_functional_equivalence_mode ? false : (params.dragen_maximum_quality_mode ? true : params.use_spanning_event_genotyping)
    def unmap_contaminant_reads_ = params.dragen_functional_equivalence_mode ? false : (params.dragen_maximum_quality_mode ? true : params.unmap_contaminant_reads)
    def perform_bqsr_ = (params.dragen_functional_equivalence_mode || params.dragen_maximum_quality_mode) ? false : params.perform_bqsr
    def use_gatk3_haplotype_caller_ = (params.dragen_functional_equivalence_mode || params.dragen_maximum_quality_mode) ? false : params.use_gatk3_haplotype_caller

    // Validate required parameters
    if (!params.input_fastq_r1) error "Missing required parameter: input_fastq_r1"
    if (!params.input_fastq_r2) error "Missing required parameter: input_fastq_r2"
    if (!params.sample_name) error "Missing required parameter: sample_name"
    if (!params.reference_fasta) error "Missing required parameter: reference_fasta"

    // Validate conflicting parameters
    if (params.dragen_functional_equivalence_mode && params.dragen_maximum_quality_mode) {
        error "Both dragen_functional_equivalence_mode and dragen_maximum_quality_mode have been set to true, however, they are mutually exclusive."
    }

    if (run_dragen_mode_variant_calling_ && use_gatk3_haplotype_caller_) {
        error "DRAGEN mode variant calling has been activated, however, the HaplotypeCaller version has been set to use GATK 3."
    }

    // Create input channels
    // Reference files channels
    reference_fasta_ch = Channel.fromPath(params.reference_fasta, checkIfExists: true)
    reference_fasta_index_ch = Channel.fromPath(params.reference_fasta_index, checkIfExists: true)
    reference_dict_ch = Channel.fromPath(params.reference_dict, checkIfExists: true)

    str_table_ch = Channel.fromPath(params.str_table_file, checkIfExists: true)

    // BWA index files
    reference_alt_ch = Channel.fromPath(params.reference_alt, checkIfExists: true)
    reference_amb_ch = Channel.fromPath(params.reference_amb, checkIfExists: true)
    reference_ann_ch = Channel.fromPath(params.reference_ann, checkIfExists: true)
    reference_bwt_ch = Channel.fromPath(params.reference_bwt, checkIfExists: true)
    reference_pac_ch = Channel.fromPath(params.reference_pac, checkIfExists: true)
    reference_sa_ch = Channel.fromPath(params.reference_sa, checkIfExists: true)

    // FASTQ files
    fastq_r1_ch = Channel.fromPath(params.input_fastq_r1, checkIfExists: true)
    fastq_r2_ch = Channel.fromPath(params.input_fastq_r2, checkIfExists: true)
    fastq_pairs_ch = fastq_r1_ch.combine(fastq_r2_ch)

    // Initialize variables to avoid undefined issues
    aligned_bam_ch = Channel.empty()
    unmapped_bam_ch = Channel.empty()

    // Alignment - conditional based on aligner parameter
    if (params.aligner == "dragmap") {
        // Validate DRAGMAP reference files
        if (!params.dragmap_reference_bin || !params.dragmap_hash_table_cfg_bin || !params.dragmap_hash_table_cmp) {
            error "DRAGMAP alignment selected but required DRAGMAP reference files are missing: dragmap_reference_bin, dragmap_hash_table_cfg_bin, dragmap_hash_table_cmp"
        }

        // Convert FASTQ to unmapped BAM
        FASTQ2UBAM(
            fastq_r1_ch,
            fastq_r2_ch,
            params.sample_name,
            params.read_group_id ?: params.sample_name,
            params.read_group_platform,
            params.read_group_pu ?: "unknown",
            params.read_group_library ?: params.sample_name,
            params.read_group_center ?: "unknown"
        )
        unmapped_bam_ch = FASTQ2UBAM.out.unmapped_bam

        // DRAGMAP reference files channels
        dragmap_reference_bin_ch = Channel.fromPath(params.dragmap_reference_bin, checkIfExists: true)
        dragmap_hash_table_cfg_bin_ch = Channel.fromPath(params.dragmap_hash_table_cfg_bin, checkIfExists: true)
        dragmap_hash_table_cmp_ch = Channel.fromPath(params.dragmap_hash_table_cmp, checkIfExists: true)

        DRAGMAP_ALIGN(
            unmapped_bam_ch,
            dragmap_reference_bin_ch,
            dragmap_hash_table_cfg_bin_ch,
            dragmap_hash_table_cmp_ch,
            reference_fasta_ch,
            reference_dict_ch,
            params.sample_name,
            params.read_group_id ?: params.sample_name,
            params.read_group_platform,
            params.read_group_pu ?: "unknown",
            params.read_group_library ?: params.sample_name,
            params.read_group_center ?: "unknown"
        )

        aligned_bam_ch = DRAGMAP_ALIGN.out.aligned_bam
    } else if (params.aligner == "bwa-mem") {
        // BWA-MEM alignment
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

        aligned_bam_ch = BWA_MEM_ALIGN.out.aligned_bam
    } else {
        error "Invalid aligner specified: '${params.aligner}'. Valid options are 'dragmap' or 'bwa-mem'"
    }

    // Mark duplicates
    MARK_DUPLICATES(
        aligned_bam_ch,
        base_file_name
    )

    // Sort BAM
    SORT_BAM(
        MARK_DUPLICATES.out.marked_bam,
        base_file_name
    )

    // Dragen dragstr model calibration (if enabled)
    if (run_dragen_mode_variant_calling_) {
            CALIBRATE_DRAGSTR_MODEL(
            SORT_BAM.out.sorted_bam,
            SORT_BAM.out.sorted_bai,
            reference_fasta_ch,
            reference_fasta_index_ch,
            reference_dict_ch,
            str_table_ch,
            params.sample_name
        )
    }

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
            reference_fasta_ch,
            reference_fasta_index_ch,
            reference_dict_ch,
            BASE_RECALIBRATOR.out.recal_table,
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
        reference_dict_ch,
        wgs_interval_ch,
        base_file_name,
        read_length
    )

    COLLECT_RAW_WGS_METRICS(
        final_bam_ch,
        final_bam_index_ch,
        reference_fasta_ch,
        reference_fasta_index_ch,
        reference_dict_ch,
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
        reference_fasta_ch,
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
        use_spanning_event_genotyping_,
        CALIBRATE_DRAGSTR_MODEL.out.dragstr_model
    )

    if (run_dragen_mode_variant_calling_) {

        DRAGEN_HARD_VARIANT_FILTRATION(
            HAPLOTYPE_CALLER.out.gvcf,
            HAPLOTYPE_CALLER.out.gvcf_index,
            reference_fasta_ch,
            reference_fasta_index_ch,
            reference_dict_ch,
            base_file_name
        )
    }

}

workflow.onComplete {

    println "Pipeline completed: ${workflow.success ? 'OK' : 'failed'}"
    if (workflow.success) {

        println """
        Results are located in: ${params.outdir}
        
        Main outputs:
        - GVCF: ${params.outdir}/variants/${params.final_gvcf_base_name ?: params.base_file_name ?: params.sample_name}.g.vcf.gz
        - CRAM: ${params.outdir}/cram/${params.base_file_name ?: params.sample_name}.cram
        - Metrics: Various QC metrics files in ${params.outdir}/qc
        """
    }
}

workflow.onError {
        println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
    }
