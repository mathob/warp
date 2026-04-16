/*
 * Quality control processes
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
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8' }"
    
    input:
    path input_bam
    path input_bai
    path reference_fasta
    path reference_fasta_index
    path reference_dict
    val base_name
    
    output:
    path "${base_name}.wgs_metrics.txt", emit: wgs_metrics
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def memory_gb = task.memory.toGiga()
    
    """
    java -Xmx${memory_gb-1}G -jar \${PICARD_JAR} CollectWgsMetrics \\
        INPUT=${input_bam} \\
        OUTPUT=${base_name}.wgs_metrics.txt \\
        REFERENCE_SEQUENCE=${reference_fasta} \\
        MINIMUM_MAPPING_QUALITY=20 \\
        MINIMUM_BASE_QUALITY=20 \\
        COVERAGE_CAP=250 \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(java -jar \${PICARD_JAR} CollectWgsMetrics --version 2>&1) | grep -o 'Version:[0-9.]*' | cut -f2 -d:)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.wgs_metrics.txt
    touch versions.yml
    """
}

/*
 * Collect raw WGS metrics using Picard
 */
process COLLECT_RAW_WGS_METRICS {
    tag "${base_name}"
    label 'process_medium'
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::picard=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8' }"
    
    input:
    path input_bam
    path input_bai
    path reference_fasta
    path reference_fasta_index
    path reference_dict
    val base_name
    
    output:
    path "${base_name}.raw_wgs_metrics.txt", emit: raw_wgs_metrics
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def memory_gb = task.memory.toGiga()
    
    """
    java -Xmx${memory_gb-1}G -jar \${PICARD_JAR} CollectRawWgsMetrics \\
        INPUT=${input_bam} \\
        OUTPUT=${base_name}.raw_wgs_metrics.txt \\
        REFERENCE_SEQUENCE=${reference_fasta} \\
        MINIMUM_MAPPING_QUALITY=20 \\
        MINIMUM_BASE_QUALITY=20 \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(java -jar \${PICARD_JAR} CollectRawWgsMetrics --version 2>&1) | grep -o 'Version:[0-9.]*' | cut -f2 -d:)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.raw_wgs_metrics.txt
    touch versions.yml
    """
}

/*
 * Check contamination using VerifyBamID2
 */
process CHECK_CONTAMINATION {
    tag "${base_name}"
    label 'process_medium'
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::verifybamid2=2.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/verify-bam-id:1.0.0-c1cba76e979904eb69c31520a0d7f5be63c72253-1626442707' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/verify-bam-id:1.0.0-c1cba76e979904eb69c31520a0d7f5be63c72253-1626442707' }"
    
    input:
    path input_bam
    path input_bai
    path contamination_sites_bed
    path contamination_sites_mu
    path contamination_sites_ud
    val base_name
    
    output:
    path "${base_name}.selfSM", emit: contamination_selfSM
    path "${base_name}.Ancestry", emit: contamination_ancestry, optional: true
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    
    """
    VerifyBamID \\
        --UDPath ${contamination_sites_ud} \\
        --MeanPath ${contamination_sites_mu} \\
        --BedPath ${contamination_sites_bed} \\
        --BamFile ${input_bam} \\
        --Output ${base_name} \\
        --DisableSanityCheck \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verifybamid2: \$(VerifyBamID --help | head -1 | sed 's/^.*VerifyBamID //; s/ .*\$//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.selfSM
    touch ${base_name}.Ancestry
    touch versions.yml
    """
}

/*
 * Collect aggregated BAM QC metrics
 */
process AGGREGATED_BAM_QC {
    tag "${base_name}"
    label 'process_low'
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::picard=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8' }"
    
    input:
    path duplicate_metrics
    path wgs_metrics
    path raw_wgs_metrics
    path contamination_selfSM
    val base_name
    
    output:
    path "${base_name}.aggregated_metrics.txt", emit: aggregated_metrics
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    # Create aggregated metrics file
    echo "# Aggregated QC Metrics for ${base_name}" > ${base_name}.aggregated_metrics.txt
    echo "" >> ${base_name}.aggregated_metrics.txt
    
    # Add duplicate metrics
    echo "## Duplicate Metrics" >> ${base_name}.aggregated_metrics.txt
    cat ${duplicate_metrics} >> ${base_name}.aggregated_metrics.txt
    echo "" >> ${base_name}.aggregated_metrics.txt
    
    # Add WGS metrics
    echo "## WGS Metrics" >> ${base_name}.aggregated_metrics.txt
    cat ${wgs_metrics} >> ${base_name}.aggregated_metrics.txt
    echo "" >> ${base_name}.aggregated_metrics.txt
    
    # Add raw WGS metrics
    echo "## Raw WGS Metrics" >> ${base_name}.aggregated_metrics.txt
    cat ${raw_wgs_metrics} >> ${base_name}.aggregated_metrics.txt
    echo "" >> ${base_name}.aggregated_metrics.txt
    
    # Add contamination check
    echo "## Contamination Check" >> ${base_name}.aggregated_metrics.txt
    cat ${contamination_selfSM} >> ${base_name}.aggregated_metrics.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version 2>&1) | sed 's/^.*bash, version //; s/ .*\$//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.aggregated_metrics.txt
    touch versions.yml
    """
}

/*
 * Validate SAM file
 */
process VALIDATE_SAM_FILE {
    tag "${base_name}"
    label 'process_medium'
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::picard=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8' }"
    
    input:
    path input_bam
    path input_bai
    path reference_fasta
    path reference_fasta_index
    path reference_dict
    val base_name
    
    output:
    path "${base_name}.validation_report.txt", emit: validation_report
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def memory_gb = task.memory.toGiga()
    
    """
    java -Xmx${memory_gb-1}G -jar \${PICARD_JAR} ValidateSamFile \\
        INPUT=${input_bam} \\
        OUTPUT=${base_name}.validation_report.txt \\
        REFERENCE_SEQUENCE=${reference_fasta} \\
        MODE=SUMMARY \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(java -jar \${PICARD_JAR} ValidateSamFile --version 2>&1) | grep -o 'Version:[0-9.]*' | cut -f2 -d:)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.validation_report.txt
    touch versions.yml
    """
}
