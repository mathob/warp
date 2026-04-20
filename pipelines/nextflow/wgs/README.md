# WholeGenomeGermlineSingleSampleFastq - Nextflow Implementation

This directory contains a Nextflow implementation of the WARP WholeGenomeGermlineSingleSampleFastq pipeline, originally written in WDL.

## Overview

This pipeline implements data pre-processing and initial variant calling (GVCF generation) according to the GATK Best Practices for germline SNP and Indel discovery in human whole-genome data.

**Original WDL Pipeline**: https://github.com/broadinstitute/warp/blob/master/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSampleFastq.wdl

## Pipeline Steps

1. **Alignment**: BWA-MEM alignment of paired-end FASTQ files to reference genome
2. **Duplicate Marking**: Mark PCR/optical duplicates using Picard MarkDuplicates
3. **Sorting**: Sort BAM files by coordinate
4. **Base Quality Score Recalibration (BQSR)**: 
   - Generate recalibration table with GATK BaseRecalibrator
   - Apply recalibration with GATK ApplyBQSR
5. **Quality Control**: Generate comprehensive QC metrics
6. **Contamination Check**: Estimate sample contamination using VerifyBamID2
7. **Format Conversion**: Convert final BAM to CRAM format
8. **Variant Calling**: Generate GVCF using GATK HaplotypeCaller

## Requirements

### Input Requirements
- Human whole-genome paired-end sequencing data in FASTQ format
- Reference genome files (Hg38 with ALT contigs recommended)
- Known variant sites for base recalibration
- Contamination reference files

### Software Requirements
- Nextflow (>=21.10.3)
- One of the following:
  - Docker
  - Singularity/Apptainer
  - Conda/Mamba
  - Manual installation of required tools

### Required Tools
- BWA (>=0.7.17)
- GATK4 (>=4.4.0.0)
- Picard (>=3.0.0)
- Samtools (>=1.17)
- VerifyBamID2 (>=2.0.1)

## Quick Start

### 1. Clone the Repository
```bash
git clone https://github.com/broadinstitute/warp.git
cd warp/pipelines/nextflow/wgs
```

### 2. Prepare Your Parameters File
Copy and edit the example parameters file:
```bash
cp params.config my_params.config
# Edit my_params.config with your file paths and settings
```

### 3. Run the Pipeline

**With Docker:**
```bash
nextflow run WholeGenomeGermlineSingleSampleFastq.nf \
    -c my_params.config \
    -profile docker
```

**With Singularity:**
```bash
nextflow run WholeGenomeGermlineSingleSampleFastq.nf \
    -c my_params.config \
    -profile singularity
```

**With Conda:**
```bash
nextflow run WholeGenomeGermlineSingleSampleFastq.nf \
    -c my_params.config \
    -profile conda
```

**On SLURM cluster:**
```bash
nextflow run WholeGenomeGermlineSingleSampleFastq.nf \
    -c my_params.config \
    -profile slurm,singularity
```

## Configuration

### Required Parameters

The following parameters must be specified in your configuration file:

```groovy
params {
    // Input files
    input_fastq_r1 = "/path/to/sample_R1.fastq.gz"
    input_fastq_r2 = "/path/to/sample_R2.fastq.gz"
    sample_name = "your_sample_name"
    
    // Reference files
    reference_fasta = "/path/to/Homo_sapiens_assembly38.fasta"
    reference_fasta_index = "/path/to/Homo_sapiens_assembly38.fasta.fai"
    reference_dict = "/path/to/Homo_sapiens_assembly38.dict"
    
    // BWA index files
    reference_alt = "/path/to/Homo_sapiens_assembly38.fasta.64.alt"
    reference_amb = "/path/to/Homo_sapiens_assembly38.fasta.64.amb"
    reference_ann = "/path/to/Homo_sapiens_assembly38.fasta.64.ann"
    reference_bwt = "/path/to/Homo_sapiens_assembly38.fasta.64.bwt"
    reference_pac = "/path/to/Homo_sapiens_assembly38.fasta.64.pac"
    reference_sa = "/path/to/Homo_sapiens_assembly38.fasta.64.sa"
    
    // Additional required files...
    // See params.config for complete list
}
```

### Optional Parameters

Many parameters have sensible defaults but can be customized:

- `provide_bam_output`: Set to `true` if you need final BAM files (default: `false`)
- `perform_bqsr`: Enable/disable base quality score recalibration (default: `true`)
- `aligner`: Alignment algorithm to use ("dragmap" or "bwa-mem", default: "dragmap")
- `dragen_functional_equivalence_mode`: DRAGEN compatibility mode (default: `false`)

## Output Files

Results are organized in the output directory (default: `results/`):

```
results/
├── cram/                           # CRAM files and indices
│   ├── {sample}.cram
│   ├── {sample}.cram.crai
│   └── {sample}.cram.md5
├── vcf/                           # Variant calls
│   ├── {sample}.g.vcf.gz
│   └── {sample}.g.vcf.gz.tbi
├── qc/                            # Quality control metrics
│   ├── {sample}.wgs_metrics.txt
│   ├── {sample}.raw_wgs_metrics.txt
│   ├── {sample}.selfSM
│   └── {sample}.contamination.txt
└── pipeline_info/                 # Pipeline execution reports
    ├── execution_report.html
    ├── execution_timeline.html
    └── execution_trace.txt
```

## Execution Profiles

Several execution profiles are available:

- `docker`: Use Docker containers
- `singularity`: Use Singularity/Apptainer containers  
- `conda`: Use Conda environments
- `slurm`: Submit jobs to SLURM scheduler
- `pbs`: Submit jobs to PBS scheduler
- `lsf`: Submit jobs to LSF scheduler
- `aws`: Run on AWS Batch
- `gcp`: Run on Google Cloud Batch

Profiles can be combined, e.g., `-profile slurm,singularity`

## Resource Requirements

### Minimum Resources
- CPU: 8 cores
- Memory: 32 GB RAM  
- Storage: ~500 GB for 30x WGS sample

### Recommended Resources
- CPU: 16+ cores
- Memory: 64+ GB RAM
- Storage: 1+ TB for multiple samples

### Execution Time
- ~24-48 hours for 30x WGS sample (depending on resources)

## Troubleshooting

### Common Issues

1. **Out of Memory Errors**
   - Increase memory allocation in nextflow.config
   - Use fewer parallel processes

2. **Missing Reference Files**
   - Ensure all reference file paths are correct
   - Check file permissions and accessibility

3. **Container Issues**
   - Verify Docker/Singularity is properly installed
   - Check container registry accessibility

### Getting Help

1. Check Nextflow documentation: https://nextflow.io/docs/latest/
2. Review the original WDL pipeline documentation
3. File issues on the WARP GitHub repository

## Differences from WDL Version

This Nextflow implementation aims to maintain functional equivalence with the original WDL pipeline while adapting to Nextflow's execution model:

### Key Differences:
- **Execution Model**: Nextflow's dataflow model vs WDL's call-based model
- **Parallelization**: Some scatter-gather operations simplified for initial version
- **Error Handling**: Adapted to Nextflow's error handling mechanisms
- **Resource Management**: Uses Nextflow's process-based resource allocation

### Maintained Features:
- Same tools and versions
- Identical parameter logic
- Same quality control metrics
- Compatible outputs

## Contributing

To contribute improvements or bug fixes:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

This pipeline is released under the same license as the original WARP pipelines. See the main repository for license details.
