/*
 * Quality control processes for WGS data
 */

/*
 * Collect WGS metrics using Picard
 */
process COLLECT_WGS_METRICS {
    tag "${base_name}"
    label 'process_medium'
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::picard=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'quay.io/biocontainers/picard:3.0.0--hdfd78af_1' }"
    
    input:
    path input_bam
    path input_bai
    path reference_fasta
    path reference_fasta_index
    path wgs_coverage_interval_list
    val base_name
    val read_length
    
    output:
    path "${base_name}.wgs_metrics.txt", emit: metrics
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectWgsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    
    """
    java -Xmx${avail_mem}M -jar \$PICARD_HOME/picard.jar CollectWgsMetrics \\
        INPUT=${input_bam} \\
        OUTPUT=${base_name}.wgs_metrics.txt \\
        REFERENCE_SEQUENCE=${reference_fasta} \\
        INTERVALS=${wgs_coverage_interval_list} \\
        MINIMUM_MAPPING_QUALITY=20 \\
        MINIMUM_BASE_QUALITY=20 \\
        COVERAGE_CAP=250 \\
        READ_LENGTH=${read_length} \\
        VALIDATION_STRINGENCY=SILENT \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(java -jar \$PICARD_HOME/picard.jar CollectWgsMetrics --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.wgs_metrics.txt
    touch versions.yml
    """
}

/*
 * Collect raw WGS metrics using Picard (less stringent thresholds)
 */
process COLLECT_RAW_WGS_METRICS {
    tag "${base_name}"
    label 'process_medium'
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::picard=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'quay.io/biocontainers/picard:3.0.0--hdfd78af_1' }"
    
    input:
    path input_bam
    path input_bai
    path reference_fasta
    path reference_fasta_index
    path wgs_coverage_interval_list
    val base_name
    val read_length
    
    output:
    path "${base_name}.raw_wgs_metrics.txt", emit: metrics
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectRawWgsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    
    """
    java -Xmx${avail_mem}M -jar \$PICARD_HOME/picard.jar CollectWgsMetrics \\
        INPUT=${input_bam} \\
        OUTPUT=${base_name}.raw_wgs_metrics.txt \\
        REFERENCE_SEQUENCE=${reference_fasta} \\
        INTERVALS=${wgs_coverage_interval_list} \\
        MINIMUM_MAPPING_QUALITY=3 \\
        MINIMUM_BASE_QUALITY=3 \\
        COVERAGE_CAP=250 \\
        READ_LENGTH=${read_length} \\
        VALIDATION_STRINGENCY=SILENT \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(java -jar \$PICARD_HOME/picard.jar CollectWgsMetrics --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.raw_wgs_metrics.txt
    touch versions.yml
    """
}

/*
 * Check contamination using VerifyBamID
 */
process CHECK_CONTAMINATION {
    tag "${base_name}"
    label 'process_medium'
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::verifybamid2=2.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verifybamid2:2.0.1--h1b792b2_0' :
        'quay.io/biocontainers/verifybamid2:2.0.1--h1b792b2_0' }"
    
    input:
    path input_bam
    path input_bai
    path contamination_sites_ud
    path contamination_sites_bed
    path contamination_sites_mu
    val base_name
    
    output:
    path "${base_name}.selfSM", emit: selfSM
    path "${base_name}.contamination.txt", emit: contamination_value
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    
    """
    # Run VerifyBamID2 for contamination estimation
    VerifyBamID \\
        --Verbose \\
        --NumPC 4 \\
        --Output ${base_name} \\
        --BamFile ${input_bam} \\
        --Reference ${contamination_sites_ud} \\
        --UDPath ${contamination_sites_ud} \\
        --MeanPath ${contamination_sites_mu} \\
        --BedPath ${contamination_sites_bed} \\
        ${args}
    
    # Extract contamination fraction from selfSM file
    if [ -f "${base_name}.selfSM" ]; then
        CONTAMINATION=\$(tail -n 1 ${base_name}.selfSM | cut -f7)
        echo \$CONTAMINATION > ${base_name}.contamination.txt
    else
        echo "0.0" > ${base_name}.contamination.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verifybamid2: \$(VerifyBamID --help 2>&1 | grep -o "VerifyBamID [0-9.]*" | sed 's/VerifyBamID //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.selfSM
    echo "0.0" > ${base_name}.contamination.txt
    touch versions.yml
    """
}

/*
 * Aggregated BAM QC using multiple Picard tools
 */
process AGGREGATED_BAM_QC {
    tag "${base_name}"
    label 'process_medium'
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::picard=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'quay.io/biocontainers/picard:3.0.0--hdfd78af_1' }"
    
    input:
    path input_bam
    path input_bai
    path reference_fasta
    path reference_fasta_index
    path reference_dict
    path haplotype_database_file
    val base_name
    val sample_name
    
    output:
    path "${base_name}.alignment_summary_metrics.txt", emit: alignment_summary_metrics
    path "${base_name}.insert_size_metrics.txt", emit: insert_size_metrics
    path "${base_name}.insert_size_histogram.pdf", emit: insert_size_histogram
    path "${base_name}.gc_bias_metrics.txt", emit: gc_bias_metrics
    path "${base_name}.gc_bias_summary_metrics.txt", emit: gc_bias_summary_metrics
    path "${base_name}.gc_bias.pdf", emit: gc_bias_pdf
    path "${base_name}.quality_distribution_metrics.txt", emit: quality_distribution_metrics
    path "${base_name}.quality_distribution.pdf", emit: quality_distribution_pdf
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectMultipleMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    
    """
    # Collect multiple QC metrics
    java -Xmx${avail_mem}M -jar \$PICARD_HOME/picard.jar CollectMultipleMetrics \\
        INPUT=${input_bam} \\
        OUTPUT=${base_name} \\
        REFERENCE_SEQUENCE=${reference_fasta} \\
        PROGRAM=CollectAlignmentSummaryMetrics \\
        PROGRAM=CollectInsertSizeMetrics \\
        PROGRAM=CollectGcBiasMetrics \\
        PROGRAM=CollectSequencingArtifactMetrics \\
        PROGRAM=QualityScoreDistribution \\
        METRIC_ACCUMULATION_LEVEL=null \\
        METRIC_ACCUMULATION_LEVEL=SAMPLE \\
        METRIC_ACCUMULATION_LEVEL=LIBRARY \\
        VALIDATION_STRINGENCY=SILENT \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(java -jar \$PICARD_HOME/picard.jar CollectMultipleMetrics --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.alignment_summary_metrics.txt
    touch ${base_name}.insert_size_metrics.txt
    touch ${base_name}.insert_size_histogram.pdf
    touch ${base_name}.gc_bias_metrics.txt
    touch ${base_name}.gc_bias_summary_metrics.txt
    touch ${base_name}.gc_bias.pdf
    touch ${base_name}.quality_distribution_metrics.txt
    touch ${base_name}.quality_distribution.pdf
    touch versions.yml
    """
}

/*
 * Validate BAM file using Picard
 */
process VALIDATE_SAM_FILE {
    tag "${base_name}"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::picard=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'quay.io/biocontainers/picard:3.0.0--hdfd78af_1' }"
    
    input:
    path input_bam
    val base_name
    
    output:
    path "${base_name}.validation_report.txt", emit: validation_report
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard ValidateSamFile] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    
    """
    java -Xmx${avail_mem}M -jar \$PICARD_HOME/picard.jar ValidateSamFile \\
        INPUT=${input_bam} \\
        OUTPUT=${base_name}.validation_report.txt \\
        REFERENCE_SEQUENCE=null \\
        MAX_OUTPUT=100 \\
        IGNORE_WARNINGS=false \\
        SKIP_MATE_VALIDATION=false \\
        IS_BISULFITE_SEQUENCED=false \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(java -jar \$PICARD_HOME/picard.jar ValidateSamFile --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.validation_report.txt
    touch versions.yml
    """
}
