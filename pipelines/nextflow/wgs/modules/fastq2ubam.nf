/*
 * FASTQ to unmapped BAM conversion processes
 */

/*
 * Convert paired FASTQ files to unmapped BAM using Picard FastqToSam
 */
process FASTQ2UBAM {
    tag "${sample_name}"
    label 'process_medium'
    publishDir "${params.outdir}/ubam", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::picard=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10' :
        'us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10' }"
    
    input:
    path fastq_r1
    path fastq_r2
    val sample_name
    val read_group_id
    val read_group_platform
    val read_group_pu
    val read_group_library
    val read_group_center
    
    output:
    path "${sample_name}.unmapped.bam", emit: unmapped_bam
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def memory_mb = task.memory ? task.memory.toMega() - 512 : 3584
    
    """
    java -Xmx${memory_mb}m -jar \$PICARD_JAR FastqToSam \\
        FASTQ=${fastq_r1} \\
        FASTQ2=${fastq_r2} \\
        OUTPUT=${sample_name}.unmapped.bam \\
        READ_GROUP_NAME=${read_group_id} \\
        SAMPLE_NAME=${sample_name} \\
        LIBRARY_NAME=${read_group_library} \\
        PLATFORM_UNIT=${read_group_pu} \\
        PLATFORM=${read_group_platform} \\
        SEQUENCING_CENTER=${read_group_center} \\
        TMP_DIR=. \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(java -jar \$PICARD_JAR FastqToSam --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d' ')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_name}.unmapped.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 2.26.10
    END_VERSIONS
    """
}