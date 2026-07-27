version 1.0

## Copyright Broad Institute, 2018
##
## Modified: alignment-only variant of WholeGenomeGermlineSingleSample.wdl.
## Runs UnmappedBamToAlignedBam only — no AggregatedBamQC, BamToCram,
## CollectWgsMetrics, CollectRawWgsMetrics, or VariantCalling. Produces the
## coordinate-sorted, duplicate-marked, base-recalibrated BAM plus the
## per-readgroup and contamination metrics that UnmappedBamToAlignedBam emits
## internally.
##
## Requirements/expectations :
## - Human whole-genome paired-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
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


  String pipeline_version = "3.3.1-alignment-only"


  input {
    SampleAndUnmappedBams sample_and_unmapped_bams
    DNASeqSingleSampleReferences references
    DragmapReference? dragmap_reference
    PapiSettings papi_settings

    Boolean dragen_functional_equivalence_mode = false
    Boolean dragen_maximum_quality_mode = false

    Boolean unmap_contaminant_reads = true
    Boolean perform_bqsr = true
    Boolean use_bwa_mem = true
    Boolean allow_empty_ref_alt = false
  }

  if (dragen_functional_equivalence_mode && dragen_maximum_quality_mode) {
    call Utilities.ErrorWithMessage as PresetArgumentsError {
      input:
        message = "Both dragen_functional_equivalence_mode and dragen_maximum_quality_mode have been set to true, however, they are mutually exclusive. You can set either of them to true, or set them both to false and adjust the arguments individually."
    }
  }

  # Set DRAGEN-related alignment arguments according to the preset modes.
  # (Variant-calling-only presets have been dropped along with the variant-calling stage.)
  Boolean unmap_contaminant_reads_ = if dragen_functional_equivalence_mode then false else (if dragen_maximum_quality_mode then true else unmap_contaminant_reads)
  Boolean perform_bqsr_ = if (dragen_functional_equivalence_mode || dragen_maximum_quality_mode) then false else perform_bqsr
  Boolean use_bwa_mem_ = if (dragen_functional_equivalence_mode || dragen_maximum_quality_mode) then false else use_bwa_mem

  # Not overridable:
  Float lod_threshold = -20.0
  String cross_check_fingerprints_by = "READGROUP"
  String recalibrated_bam_basename = sample_and_unmapped_bams.base_file_name + ".aligned.duplicates_marked.recalibrated"

  call ToBam.UnmappedBamToAlignedBam {
    input:
      sample_and_unmapped_bams    = sample_and_unmapped_bams,
      references                  = references,
      dragmap_reference           = dragmap_reference,
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

  # Outputs that will be retained when execution is complete.
  # Only the artefacts produced by UnmappedBamToAlignedBam are exposed —
  # downstream QC/CRAM/GVCF outputs from the parent pipeline are intentionally omitted.
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
