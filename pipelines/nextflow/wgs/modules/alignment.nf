/*
 * Alignment processes for whole genome sequencing
 */

/*
 * BWA MEM alignment of paired-end FASTQ files
 */
process BWA_MEM_ALIGN {
    tag "${sample_name}"
    label 'process_high'
    
    conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.17 bioconda::picard=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438' }"
    
    input:
    tuple path(fastq_r1), path(fastq_r2)
    path reference_fasta
    path reference_fasta_index
    path reference_alt
    path reference_amb
    path reference_ann
    path reference_bwt
    path reference_pac
    path reference_sa
    val sample_name
    val read_group_id
    val read_group_platform
    val read_group_pu
    val read_group_library
    val read_group_center
    
    output:
    path "${sample_name}.aligned.bam", emit: aligned_bam
    path "${sample_name}.aligned.bam.bai", emit: aligned_bai
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def read_group = "@RG\\tID:${read_group_id}\\tSM:${sample_name}\\tPL:${read_group_platform}\\tPU:${read_group_pu}\\tLB:${read_group_library}\\tCN:${read_group_center}"
    
    """
    # Check if reference files are properly linked
    if [ ! -f ${reference_fasta} ]; then
        echo "Error: Reference FASTA not found"
        exit 1
    fi

    if [ -f /usr/gitc/bwa ]; then
        BWA=/usr/gitc/bwa
    else
        BWA=bwa
    fi
    
    # Perform BWA MEM alignment
    "\$BWA" mem \\
        -t ${task.cpus} \\
        -R "${read_group}" \\
        ${args} \\
        ${reference_fasta} \\
        ${fastq_r1} \\
        ${fastq_r2} | \\
    samtools view -@ ${task.cpus} -Sb - > temp.bam
    
    # Sort the BAM file and create index
    samtools sort -@ ${task.cpus} -o ${sample_name}.aligned.bam temp.bam
    samtools index ${sample_name}.aligned.bam
    
    # Clean up temporary file
    rm temp.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_name}.aligned.bam
    touch ${sample_name}.aligned.bam.bai
    touch versions.yml
    """
}

/*
 * DRAGMAP alignment (alternative to BWA MEM)
 */
process DRAGMAP_ALIGN {
    tag "${sample_name}"
    label 'process_high'
    
    conda (params.enable_conda ? "bioconda::dragmap=1.2.1 bioconda::samtools=1.17" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/dragmap:1.2.1' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/dragmap:1.2.1' }"
    
    input:
    tuple path(fastq_r1), path(fastq_r2)
    path dragmap_reference_bin
    path dragmap_hash_table_cfg_bin
    path dragmap_hash_table_cmp
    val sample_name
    val read_group_id
    val read_group_platform
    val read_group_pu
    val read_group_library
    val read_group_center
    
    output:
    path "${sample_name}.aligned.bam", emit: aligned_bam
    path "${sample_name}.aligned.bam.bai", emit: aligned_bai
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def read_group = "@RG\\tID:${read_group_id}\\tSM:${sample_name}\\tPL:${read_group_platform}\\tPU:${read_group_pu}\\tLB:${read_group_library}\\tCN:${read_group_center}"
    
    """
    # DRAGMAP alignment
    dragen-os \\
        --fastq-file1 ${fastq_r1} \\
        --fastq-file2 ${fastq_r2} \\
        --ref-dir . \\
        --RGID ${read_group_id} \\
        --RGSM ${sample_name} \\
        --num-threads ${task.cpus} \\
        ${args} \\
        --output-file-prefix ${sample_name}.aligned \\
    
    # Index the BAM file
    samtools index ${sample_name}.aligned.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragmap: \$(dragen-os --version 2>&1 | grep -o "dragen-os [0-9.]*" | sed 's/dragen-os //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_name}.aligned.bam
    touch ${sample_name}.aligned.bam.bai
    touch versions.yml
    """
}
