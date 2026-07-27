version 1.0

## Copyright Broad Institute, 2018
##
## Modified: alignment-only, Terra-flattened variant of WholeGenomeGermlineSingleSample.wdl.
## Runs UnmappedBamToAlignedBam only — no AggregatedBamQC, BamToCram,
## CollectWgsMetrics, CollectRawWgsMetrics, or VariantCalling.
##
## Difference from WholeGenomeGermlineSingleSampleAlignmentOnly.wdl:
## every top-level workflow input is a scalar, array, or optional scalar.
## No struct-typed workflow inputs, so Terra's Inputs configuration UI displays
## one row per field. The structs UnmappedBamToAlignedBam actually needs
## (SampleAndUnmappedBams, DNASeqSingleSampleReferences, ReferenceFasta,
## PapiSettings, DragmapReference?) are constructed inside the workflow body.
##
## Requirements/expectations :
## - Human whole-genome paired-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (default ".unmapped.bam")
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
## - Reference genome must be Hg38 with ALT contigs
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "../../../../../../tasks/broad/UnmappedBamToAlignedBam.wdl" as ToBam
import "../../../../../../tasks/broad/Utilities.wdl" as Utilities
import "../../../../../../structs/dna_seq/DNASeqStructs.wdl"

# WORKFLOW DEFINITION
workflow WholeGenomeGermlineSingleSampleAlignmentOnly {


  String pipeline_version = "3.3.1-alignment-only-flat"


  input {

    # ── Sample identity ─────────────────────────────────────────────────────
    String sample_name
    String base_file_name
    Array[File] flowcell_unmapped_bams
    String unmapped_bam_suffix = ".unmapped.bam"
    String? final_gvcf_base_name

    # ── Reference FASTA (+ BWA index files, sharing the fasta basename) ─────
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File? ref_str

    # ── Other reference resources ──────────────────────────────────────────
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    File calling_interval_list
    File evaluation_interval_list
    File haplotype_database_file
    File dbsnp_vcf
    File dbsnp_vcf_index
    Array[File] known_indels_sites_vcfs
    Array[File] known_indels_sites_indices

    # ── Optional DRAGMAP hash-table files (only needed if use_bwa_mem = false) ──
    File? dragmap_reference_bin
    File? dragmap_hash_table_cfg_bin
    File? dragmap_hash_table_cmp

    # ── Runtime tuning (PapiSettings) ──────────────────────────────────────
    Int preemptible_tries = 3
    Int agg_preemptible_tries = 3

    # ── Behaviour flags ────────────────────────────────────────────────────
    Boolean dragen_functional_equivalence_mode = false
    Boolean dragen_maximum_quality_mode = false
    Boolean unmap_contaminant_reads = true
    Boolean perform_bqsr = true
    Boolean use_bwa_mem = true
    Boolean allow_empty_ref_alt = false
  }

  # ── Preset mutual-exclusion check ─────────────────────────────────────────
  if (dragen_functional_equivalence_mode && dragen_maximum_quality_mode) {
    call Utilities.ErrorWithMessage as PresetArgumentsError {
      input:
        message = "Both dragen_functional_equivalence_mode and dragen_maximum_quality_mode have been set to true, however, they are mutually exclusive. You can set either of them to true, or set them both to false and adjust the arguments individually."
    }
  }

  # ── Resolve DRAGEN-related alignment presets ──────────────────────────────
  Boolean unmap_contaminant_reads_ = if dragen_functional_equivalence_mode then false else (if dragen_maximum_quality_mode then true else unmap_contaminant_reads)
  Boolean perform_bqsr_ = if (dragen_functional_equivalence_mode || dragen_maximum_quality_mode) then false else perform_bqsr
  Boolean use_bwa_mem_ = if (dragen_functional_equivalence_mode || dragen_maximum_quality_mode) then false else use_bwa_mem

  # ── Not overridable ───────────────────────────────────────────────────────
  Float lod_threshold = -20.0
  String cross_check_fingerprints_by = "READGROUP"
  String recalibrated_bam_basename = base_file_name + ".aligned.duplicates_marked.recalibrated"

  # ── Assemble the structs UnmappedBamToAlignedBam needs ───────────────────
  ReferenceFasta reference_fasta = object {
    ref_dict: ref_dict,
    ref_fasta: ref_fasta,
    ref_fasta_index: ref_fasta_index,
    ref_alt: ref_alt,
    ref_sa: ref_sa,
    ref_amb: ref_amb,
    ref_bwt: ref_bwt,
    ref_ann: ref_ann,
    ref_pac: ref_pac,
    ref_str: ref_str
  }

  DNASeqSingleSampleReferences references = object {
    contamination_sites_ud: contamination_sites_ud,
    contamination_sites_bed: contamination_sites_bed,
    contamination_sites_mu: contamination_sites_mu,
    calling_interval_list: calling_interval_list,
    reference_fasta: reference_fasta,
    known_indels_sites_vcfs: known_indels_sites_vcfs,
    known_indels_sites_indices: known_indels_sites_indices,
    dbsnp_vcf: dbsnp_vcf,
    dbsnp_vcf_index: dbsnp_vcf_index,
    evaluation_interval_list: evaluation_interval_list,
    haplotype_database_file: haplotype_database_file
  }

  SampleAndUnmappedBams sample_and_unmapped_bams = object {
    base_file_name: base_file_name,
    final_gvcf_base_name: final_gvcf_base_name,
    flowcell_unmapped_bams: flowcell_unmapped_bams,
    sample_name: sample_name,
    unmapped_bam_suffix: unmapped_bam_suffix
  }

  PapiSettings papi_settings = object {
    preemptible_tries: preemptible_tries,
    agg_preemptible_tries: agg_preemptible_tries
  }

  # Optional DRAGMAP reference — construct only if the .bin file was supplied.
  if (defined(dragmap_reference_bin)) {
    DragmapReference dragmap_reference_constructed = object {
      reference_bin: select_first([dragmap_reference_bin]),
      hash_table_cfg_bin: select_first([dragmap_hash_table_cfg_bin]),
      hash_table_cmp: select_first([dragmap_hash_table_cmp])
    }
  }

  # ── Only call: the alignment subworkflow ─────────────────────────────────
  call ToBam.UnmappedBamToAlignedBam {
    input:
      sample_and_unmapped_bams    = sample_and_unmapped_bams,
      references                  = references,
      dragmap_reference           = dragmap_reference_constructed,
      papi_settings               = papi_settings,

      contamination_sites_ud = references.contamination_sites_ud,
      contamination_sites_bed = references.contamination_sites_bed,
      contamination_sites_mu = references.contamination_sites_mu,

      cross_check_fingerprints_by = cross_check_fingerprints_by,
      haplotype_database_file     = references.haplotype_database_file,
      lod_threshold               = lod_threshold,
      recalibrated_bam_basename   = recalibrated_bam_basename,
      perform_bqsr                = perform_bqsr_,
      use_bwa_mem                 = use_bwa_mem_,
      unmap_contaminant_reads     = unmap_contaminant_reads_,
      allow_empty_ref_alt         = allow_empty_ref_alt
  }

  # ── Outputs — only what UnmappedBamToAlignedBam produces ─────────────────
  output {
    File output_bam = UnmappedBamToAlignedBam.output_bam
    File output_bam_index = UnmappedBamToAlignedBam.output_bam_index

    File duplicate_metrics = UnmappedBamToAlignedBam.duplicate_metrics
    File? output_bqsr_reports = UnmappedBamToAlignedBam.output_bqsr_reports

    Array[File] quality_yield_metrics = UnmappedBamToAlignedBam.quality_yield_metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = UnmappedBamToAlignedBam.unsorted_read_group_base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = UnmappedBamToAlignedBam.unsorted_read_group_base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = UnmappedBamToAlignedBam.unsorted_read_group_insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = UnmappedBamToAlignedBam.unsorted_read_group_insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = UnmappedBamToAlignedBam.unsorted_read_group_quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = UnmappedBamToAlignedBam.unsorted_read_group_quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = UnmappedBamToAlignedBam.unsorted_read_group_quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = UnmappedBamToAlignedBam.unsorted_read_group_quality_distribution_metrics

    File? cross_check_fingerprints_metrics = UnmappedBamToAlignedBam.cross_check_fingerprints_metrics

    File selfSM = UnmappedBamToAlignedBam.selfSM
    Float contamination = UnmappedBamToAlignedBam.contamination
  }
  meta {
    allowNestedInputs: true
  }
}
