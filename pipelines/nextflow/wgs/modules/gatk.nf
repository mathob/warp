/*
 * GATK processes for base recalibration and variant calling
 */

/*
 * GATK BaseRecalibrator - Generate base recalibration table
 */
process BASE_RECALIBRATOR {
    tag "${base_name}"
    label 'process_medium'
    
    conda (params.enable_conda ? "bioconda::gatk4=4.4.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.1.8.0' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.1.8.0' }"
    
    input:
    path input_bam
    path reference_fasta
    path reference_fasta_index
    path reference_dict
    path known_sites
    path known_sites_indices
    val base_name
    
    output:
    path "${base_name}.recal_data.table", emit: recal_table
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def known_sites_args = known_sites.collect { "--known-sites ${it}" }.join(' ')
    def memory_gb = task.memory.toGiga()
    
    """
    gatk --java-options "-Xmx${memory_gb-1}G" BaseRecalibrator \\
        --input ${input_bam} \\
        --reference ${reference_fasta} \\
        ${known_sites_args} \\
        --output ${base_name}.recal_data.table \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.recal_data.table
    touch versions.yml
    """
}

/*
 * GATK ApplyBQSR - Apply base quality score recalibration
 */
process APPLY_BQSR {
    tag "${base_name}"
    label 'process_medium'
    
    conda (params.enable_conda ? "bioconda::gatk4=4.4.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.1.8.0' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.1.8.0' }"
    
    input:
    path input_bam
    path reference_fasta
    path reference_fasta_index
    path reference_dict
    path recal_table
    val base_name
    
    output:
    path "${base_name}.recalibrated.bam", emit: recalibrated_bam
    path "${base_name}.recalibrated.bam.bai", emit: recalibrated_bai
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def memory_gb = task.memory.toGiga()
    
    """
    gatk --java-options "-Xmx${memory_gb-1}G" ApplyBQSR \\
        --input ${input_bam} \\
        --reference ${reference_fasta} \\
        --bqsr-recal-file ${recal_table} \\
        --output ${base_name}.recalibrated.bam \\
        --create-output-bam-index true \\
        ${args}
    
    mv ${base_name}.recalibrated.bai ${base_name}.recalibrated.bam.bai
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.recalibrated.bam
    touch ${base_name}.recalibrated.bam.bai
    touch versions.yml
    """
}

/*
 * GATK HaplotypeCaller - Variant calling
 */
process HAPLOTYPE_CALLER {
    tag "${base_name}"
    label 'process_high'
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::gatk4=4.4.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.2.2.0' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.2.2.0' }"
    
    input:
    path input_bam
    path input_bai
    path reference_fasta
    path reference_fasta_index
    path reference_dict
    path calling_interval_list
    path dbsnp_vcf
    path dbsnp_vcf_index
    val contamination_value
    val base_name
    val use_gatk3_haplotype_caller
    val run_dragen_mode_variant_calling
    val use_spanning_event_genotyping
    path dragstr_model
    
    output:
    path "${base_name}.g.vcf.gz", emit: gvcf
    path "${base_name}.g.vcf.gz.tbi", emit: gvcf_index
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def memory_gb = task.memory.toGiga()
    def intervals_arg = calling_interval_list ? "--intervals ${calling_interval_list}" : ""
    def dbsnp_arg = dbsnp_vcf ? "--dbsnp ${dbsnp_vcf}" : ""
    def contamination_arg = contamination_value ? "--contamination-fraction-to-filter ${contamination_value}" : ""
    def dragen_mode_arg = run_dragen_mode_variant_calling ? "--dragen-mode" : ""
    def dragstr_model_arg = dragstr_model ? "--dragstr-params-path ${dragstr_model}" : ""
    def spanning_event_genotyping_arg = use_spanning_event_genotyping ? "" : "--disable-spanning-event-genotyping"
    
    """
    gatk --java-options "-Xmx${memory_gb-1}G" HaplotypeCaller \\
        --reference ${reference_fasta} \\
        --input ${input_bam} \\
        ${intervals_arg} \\
        --output ${base_name}.g.vcf.gz \\
        ${contamination_arg} \\
        -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \\
        ${dragen_mode_arg} \\
        ${spanning_event_genotyping_arg} \\
        -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
        ${dragstr_model_arg} \\
        --emit-ref-confidence GVCF \\
        ${dbsnp_arg} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.g.vcf.gz
    touch ${base_name}.g.vcf.gz.tbi
    touch versions.yml
    """
}

/*
 * GATK VariantFiltration - Filter variants
 */
process VARIANT_FILTRATION {
    tag "${base_name}"
    label 'process_medium'
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::gatk4=4.4.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.1.8.0' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.1.8.0' }"
    
    input:
    path input_vcf
    path input_vcf_index
    path reference_fasta
    path reference_fasta_index
    path reference_dict
    val base_name
    
    output:
    path "${base_name}.filtered.vcf.gz", emit: filtered_vcf
    path "${base_name}.filtered.vcf.gz.tbi", emit: filtered_vcf_index
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def memory_gb = task.memory.toGiga()
    
    """
    gatk --java-options "-Xmx${memory_gb-1}G" VariantFiltration \\
        --variant ${input_vcf} \\
        --reference ${reference_fasta} \\
        --output ${base_name}.filtered.vcf.gz \\
        --filter-expression "QD < 2.0" --filter-name "QD2" \\
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \\
        --filter-expression "SOR > 3.0" --filter-name "SOR3" \\
        --filter-expression "FS > 60.0" --filter-name "FS60" \\
        --filter-expression "MQ < 40.0" --filter-name "MQ40" \\
        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.filtered.vcf.gz
    touch ${base_name}.filtered.vcf.gz.tbi
    touch versions.yml
    """
}
