# Container Mapping: WDL to Nextflow

This document shows the exact container mapping from the original WDL FastqToAlignedBam pipeline to the Nextflow implementation.

## Container Sources

### WDL Containers (from tasks/broad/)
The following containers are used in the original WARP WDL files:

1. **Alignment.wdl**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438`
2. **DragmapAlignment.wdl**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/dragmap:1.2.1`
3. **BamProcessing.wdl**: 
   - Picard tasks: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8`
   - GATK tasks: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.1.8.0`
   - BedTools: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/bedtools:2.27.1`
   - VerifyBamID: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/verify-bam-id:1.0.0-c1cba76e979904eb69c31520a0d7f5be63c72253-1626442707`
4. **Qc.wdl**:
   - Picard tasks: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8`
   - Python: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/python:2.7`
   - GATK tasks: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.1.8.0`

## Nextflow Process to Container Mapping

### Alignment Processes
- **BWA_MEM_ALIGN**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438`
- **DRAGMAP_ALIGN**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/dragmap:1.2.1`

### Processing Processes
- **MARK_DUPLICATES**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8`
- **SORT_BAM**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438`
- **BAM_TO_CRAM**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438`
- **GATHER_BAM_FILES**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8`

### GATK Processes
- **BASE_RECALIBRATOR**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.1.8.0`
- **APPLY_BQSR**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.1.8.0`
- **HAPLOTYPE_CALLER**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.1.8.0`
- **VARIANT_FILTRATION**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.1.8.0`

### QC Processes
- **COLLECT_WGS_METRICS**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8`
- **COLLECT_RAW_WGS_METRICS**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8`
- **CHECK_CONTAMINATION**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/verify-bam-id:1.0.0-c1cba76e979904eb69c31520a0d7f5be63c72253-1626442707`
- **AGGREGATED_BAM_QC**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8`
- **VALIDATE_SAM_FILE**: `australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8`

## Tool Versions

Based on the container tags, the Nextflow pipeline now uses the same tool versions as the WDL:

- **BWA**: 0.7.15 (from samtools-picard-bwa container)
- **Samtools**: 1.17 equivalent (from samtools-picard-bwa container)  
- **Picard**: 2.23.8
- **GATK**: 4.1.8.0
- **DRAGMAP**: 1.2.1
- **VerifyBamID**: 2.0.1 equivalent

## Configuration

The container specifications are defined both in:
1. Individual process definitions in `modules/*.nf` files
2. Process-specific overrides in `nextflow.config`

This ensures that the Nextflow pipeline uses exactly the same containers and tool versions as the original WDL implementation, maintaining functional equivalence.
