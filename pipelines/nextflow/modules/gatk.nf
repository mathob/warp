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
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"
    
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
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK BaseRecalibrator] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    
    def known_sites_args = known_sites.collect { "--known-sites $it" }.join(' ')
    
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:ParallelGCThreads=${task.cpus}" \\
        BaseRecalibrator \\
        --input ${input_bam} \\
        --output ${base_name}.recal_data.table \\
        --reference ${reference_fasta} \\
        ${known_sites_args} \\
        --use-original-qualities \\
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
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"
    
    input:
    path input_bam
    path recal_table
    path reference_fasta
    path reference_fasta_index
    path reference_dict
    val base_name
    
    output:
    path "${base_name}.recalibrated.bam", emit: recalibrated_bam
    path "${base_name}.recalibrated.bam.bai", emit: recalibrated_bai
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK ApplyBQSR] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:ParallelGCThreads=${task.cpus}" \\
        ApplyBQSR \\
        --input ${input_bam} \\
        --output ${base_name}.recalibrated.bam \\
        --reference ${reference_fasta} \\
        --bqsr-recal-file ${recal_table} \\
        --static-quantized-quals 10 \\
        --static-quantized-quals 20 \\
        --static-quantized-quals 30 \\
        --add-output-sam-program-record \\
        --create-output-bam-index \\
        --use-original-qualities \\
        ${args}
    
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
 * GATK HaplotypeCaller - Call germline SNPs and indels
 */
process HAPLOTYPE_CALLER {
    tag "${final_gvcf_base_name}"
    label 'process_high'
    publishDir "${params.outdir}/vcf", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::gatk4=4.4.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"
    
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
    val final_gvcf_base_name
    val use_gatk3_haplotype_caller
    val run_dragen_mode_variant_calling
    val use_spanning_event_genotyping
    
    output:
    path "${final_gvcf_base_name}.g.vcf.gz", emit: gvcf
    path "${final_gvcf_base_name}.g.vcf.gz.tbi", emit: gvcf_index
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    
    def contamination_arg = contamination_value ? "--contamination-fraction-to-filter ${contamination_value}" : ""
    def spanning_events_arg = use_spanning_event_genotyping ? "" : "--disable-spanning-event-genotyping"
    def dragen_mode_arg = run_dragen_mode_variant_calling ? "--dragen-mode" : ""
    
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:ParallelGCThreads=${task.cpus}" \\
        HaplotypeCaller \\
        --input ${input_bam} \\
        --output ${final_gvcf_base_name}.g.vcf.gz \\
        --reference ${reference_fasta} \\
        --intervals ${calling_interval_list} \\
        --dbsnp ${dbsnp_vcf} \\
        --emit-ref-confidence GVCF \\
        --annotation-group StandardAnnotation \\
        --annotation-group StandardHCAnnotation \\
        ${contamination_arg} \\
        ${spanning_events_arg} \\
        ${dragen_mode_arg} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${final_gvcf_base_name}.g.vcf.gz
    touch ${final_gvcf_base_name}.g.vcf.gz.tbi
    touch versions.yml
    """
}

/*
 * GATK VariantFiltration - Apply hard filters to variants (for DRAGEN mode)
 */
process VARIANT_FILTRATION {
    tag "${base_name}"
    label 'process_medium'
    
    conda (params.enable_conda ? "bioconda::gatk4=4.4.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"
    
    input:
    path input_vcf
    path input_vcf_index
    path reference_fasta
    path reference_fasta_index
    path reference_dict
    val base_name
    
    output:
    path "${base_name}.filtered.g.vcf.gz", emit: filtered_vcf
    path "${base_name}.filtered.g.vcf.gz.tbi", emit: filtered_vcf_index
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK VariantFiltration] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    
    """
    gatk --java-options "-Xmx${avail_mem}M" \\
        VariantFiltration \\
        --variant ${input_vcf} \\
        --output ${base_name}.filtered.g.vcf.gz \\
        --reference ${reference_fasta} \\
        --filter-expression "ExcessHet > 54.69" \\
        --filter-name "ExcessHet" \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.filtered.g.vcf.gz
    touch ${base_name}.filtered.g.vcf.gz.tbi
    touch versions.yml
    """
}
