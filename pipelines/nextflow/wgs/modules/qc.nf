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
    path wgs_coverage_interval_list
    val base_name
    val read_length
    
    output:
    path "${base_name}.wgs_metrics.txt", emit: wgs_metrics
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def memory_gb = task.memory.toGiga()
    def intervals_arg = wgs_coverage_interval_list ? "INTERVALS=${wgs_coverage_interval_list}" : ""
    
    """
    java -Xmx${memory_gb-1}G -jar /usr/picard/picard.jar CollectWgsMetrics \\
        INPUT=${input_bam} \\
        OUTPUT=${base_name}.wgs_metrics.txt \\
        REFERENCE_SEQUENCE=${reference_fasta} \\
        MINIMUM_MAPPING_QUALITY=20 \\
        MINIMUM_BASE_QUALITY=20 \\
        COVERAGE_CAP=250 \\
        READ_LENGTH=${read_length} \\
        ${intervals_arg} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(java -jar /usr/picard/picard.jar CollectWgsMetrics --version 2>&1) | grep -o 'Version:[0-9.]*' | cut -f2 -d:)
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
    path wgs_coverage_interval_list
    val base_name
    val read_length
    
    output:
    path "${base_name}.raw_wgs_metrics.txt", emit: raw_wgs_metrics
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def memory_gb = task.memory.toGiga()
    def intervals_arg = wgs_coverage_interval_list ? "INTERVALS=${wgs_coverage_interval_list}" : ""
    
    """
    java -Xmx${memory_gb-1}G -jar /usr/picard/picard.jar CollectRawWgsMetrics \\
        INPUT=${input_bam} \\
        OUTPUT=${base_name}.raw_wgs_metrics.txt \\
        REFERENCE_SEQUENCE=${reference_fasta} \\
        MINIMUM_MAPPING_QUALITY=20 \\
        MINIMUM_BASE_QUALITY=20 \\
        READ_LENGTH=${read_length} \\
        ${intervals_arg} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(java -jar /usr/picard/picard.jar CollectRawWgsMetrics --version 2>&1) | grep -o 'Version:[0-9.]*' | cut -f2 -d:)
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
    path contamination_sites_ud
    path contamination_sites_bed
    path contamination_sites_mu
    val base_name
    
    output:
    path "${base_name}.selfSM", emit: selfSM
    path "${base_name}.Ancestry", emit: contamination_ancestry, optional: true
    path "${base_name}.contamination_value.txt", emit: contamination_value
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    
    """
    if [ -f /usr/gitc/VerifyBamID ]; then
        VerifyBamID=/usr/gitc/VerifyBamID
    else
        VerifyBamID=VerifyBamID
    fi
    

    "\$VerifyBamID" \\
        --UDPath ${contamination_sites_ud} \\
        --MeanPath ${contamination_sites_mu} \\
        --BedPath ${contamination_sites_bed} \\
        --BamFile ${input_bam} \\
        --Output ${base_name} \\
        --DisableSanityCheck \\
        ${args}
    
    # Extract contamination value from .selfSM file
    if [ -f ${base_name}.selfSM ]; then
        awk 'NR==2 {print \$7}' ${base_name}.selfSM > ${base_name}.contamination_value.txt
    else
        echo "0.0" > ${base_name}.contamination_value.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verifybamid2: \$(VerifyBamID --help | head -1 | sed 's/^.*VerifyBamID //; s/ .*\$//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.selfSM
    touch ${base_name}.Ancestry
    touch ${base_name}.contamination_value.txt
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
    path input_bam
    path input_bai
    path reference_fasta
    path reference_fasta_index
    path reference_dict
    path haplotype_database_file
    val base_name
    val sample_name
    
    output:
    path "${base_name}.alignment_summary_metrics", emit: alignment_summary_metrics
    path "${base_name}.insert_size_metrics", emit: insert_size_metrics
    path "${base_name}.gc_bias_metrics", emit: gc_bias_metrics
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def memory_gb = task.memory.toGiga()
    
    """
    # Collect alignment summary metrics
    java -Xmx${memory_gb-1}G -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics \\
        INPUT=${input_bam} \\
        OUTPUT=${base_name}.alignment_summary_metrics \\
        REFERENCE_SEQUENCE=${reference_fasta} \\
        ${args}
    
    # Collect insert size metrics
    java -Xmx${memory_gb-1}G -jar /usr/picard/picard.jar CollectInsertSizeMetrics \\
        INPUT=${input_bam} \\
        OUTPUT=${base_name}.insert_size_metrics \\
        HISTOGRAM_FILE=${base_name}.insert_size_histogram.pdf \\
        ${args}
    
    # Collect GC bias metrics
    java -Xmx${memory_gb-1}G -jar /usr/picard/picard.jar CollectGcBiasMetrics \\
        INPUT=${input_bam} \\
        OUTPUT=${base_name}.gc_bias_metrics \\
        CHART_OUTPUT=${base_name}.gc_bias_metrics.pdf \\
        SUMMARY_OUTPUT=${base_name}.gc_bias_summary_metrics \\
        REFERENCE_SEQUENCE=${reference_fasta} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(java -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics --version 2>&1) | grep -o 'Version:[0-9.]*' | cut -f2 -d:)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.alignment_summary_metrics
    touch ${base_name}.insert_size_metrics
    touch ${base_name}.gc_bias_metrics
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
    java -Xmx${memory_gb-1}G -jar /usr/picard/picard.jar ValidateSamFile \\
        INPUT=${input_bam} \\
        OUTPUT=${base_name}.validation_report.txt \\
        REFERENCE_SEQUENCE=${reference_fasta} \\
        MODE=SUMMARY \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(java -jar /usr/picard/picard.jar ValidateSamFile --version 2>&1) | grep -o 'Version:[0-9.]*' | cut -f2 -d:)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.validation_report.txt
    touch versions.yml
    """
}
