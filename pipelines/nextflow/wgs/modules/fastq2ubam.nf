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
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8' }"
    
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
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard FastqToSam] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    
    """
    java -Xmx${avail_mem}M -jar /usr/picard/picard.jar  FastqToSam \\
        FASTQ=${fastq_r1} \\
        FASTQ2=${fastq_r2} \\
        OUTPUT=${sample_name}.unmapped.bam \\
        READ_GROUP_NAME=${read_group_id} \\
        SAMPLE_NAME=${sample_name} \\
        LIBRARY_NAME=${read_group_library} \\
        PLATFORM_UNIT=${read_group_pu} \\
        PLATFORM=${read_group_platform} \\
        SEQUENCING_CENTER=${read_group_center} \\
        SORT_ORDER=queryname \\
        TMP_DIR=. \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(java -jar /usr/picard/picard.jar FastqToSam --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d' ')
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