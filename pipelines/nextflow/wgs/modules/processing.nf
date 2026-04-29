/*
 * BAM processing tasks (sorting, duplicate marking, format conversion)
 */

/*
 * Mark duplicates in BAM file using Picard
 */
process MARK_DUPLICATES {
    tag "${base_name}"
    label 'process_medium'
    
    conda (params.enable_conda ? "bioconda::picard=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8' }"
    
    input:
    path input_bam
    val base_name
    
    output:
    path "${base_name}.duplicates_marked.bam", emit: marked_bam
    path "${base_name}.duplicate_metrics.txt", emit: duplicate_metrics
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    
    """
    java -Xmx${avail_mem}M -jar /usr/picard/picard.jar MarkDuplicates \\
        INPUT=${input_bam} \\
        OUTPUT=${base_name}.duplicates_marked.bam \\
        METRICS_FILE=${base_name}.duplicate_metrics.txt \\
        VALIDATION_STRINGENCY=SILENT \\
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \\
        ASSUME_SORT_ORDER=queryname \\
        CLEAR_DT=false \\
        ADD_PG_TAG_TO_READS=false \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(java -jar /usr/picard/picard.jar MarkDuplicates --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.duplicates_marked.bam
    touch ${base_name}.duplicates_marked.bam.bai
    touch ${base_name}.duplicate_metrics.txt
    touch versions.yml
    """
}

/*
 * Sort BAM file by coordinate
 */
process SORT_BAM {
    tag "${base_name}"
    label 'process_medium'
    
    conda (params.enable_conda ? "bioconda::samtools=1.17" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438' }"
    
    input:
    path input_bam
    val base_name
    
    output:
    path "${base_name}.sorted.bam", emit: sorted_bam
    path "${base_name}.sorted.bam.bai", emit: sorted_bai
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -m ${task.memory.toGiga()}G \\
        ${args} \\
        -o ${base_name}.sorted.bam \\
        ${input_bam}
    
    samtools index ${base_name}.sorted.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.sorted.bam
    touch ${base_name}.sorted.bam.bai
    touch versions.yml
    """
}

/*
 * Convert BAM to CRAM format
 */
process BAM_TO_CRAM {
    tag "${base_name}"
    label 'process_medium'
    publishDir "${params.outdir}/cram", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::samtools=1.17" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/samtools-picard-bwa:1.0.0-0.7.15-2.23.8-1626449438' }"
    
    input:
    path input_bam
    path input_bai
    path reference_fasta
    path reference_fasta_index
    path reference_dict
    path duplicate_metrics
    val base_name
    
    output:
    path "${base_name}.cram", emit: cram
    path "${base_name}.cram.crai", emit: cram_index
    path "${base_name}.cram.md5", emit: cram_md5
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    
    """
    # Convert BAM to CRAM
    samtools view \\
        -@ ${task.cpus} \\
        -C \\
        -T ${reference_fasta} \\
        ${args} \\
        -o ${base_name}.cram \\
        ${input_bam}
    
    # Index CRAM file
    samtools index ${base_name}.cram
    
    # Generate MD5 checksum
    md5sum ${base_name}.cram > ${base_name}.cram.md5
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.cram
    touch ${base_name}.cram.crai
    touch ${base_name}.cram.md5
    touch versions.yml
    """
}

/*
 * Split FASTQ files into chunks for parallel alignment
 */
process SPLIT_FASTQ {

    input:
    path fastq_r1
    path fastq_r2
    val  num_chunks

    output:
    path "chunk_*_R1.fastq.gz", emit: fastq_r1_chunks
    path "chunk_*_R2.fastq.gz", emit: fastq_r2_chunks

    script:
    """
    total_reads=\$(zcat ${fastq_r1} | wc -l)
    total_reads=\$((total_reads / 4))
    reads_per_chunk=\$(( (total_reads + ${num_chunks} - 1) / ${num_chunks} ))
    lines_per_chunk=\$((reads_per_chunk * 4))

    zcat ${fastq_r1} | split -l \$lines_per_chunk --numeric-suffixes=1 --suffix-length=4 - chunk_r1_
    zcat ${fastq_r2} | split -l \$lines_per_chunk --numeric-suffixes=1 --suffix-length=4 - chunk_r2_

    for f in chunk_r1_*; do
        idx=\${f##*_}
        gzip -c \$f > chunk_\${idx}_R1.fastq.gz
        gzip -c chunk_r2_\${idx} > chunk_\${idx}_R2.fastq.gz
    done
    """

    stub:
    """
    touch chunk_0001_R1.fastq.gz chunk_0001_R2.fastq.gz
    touch chunk_0002_R1.fastq.gz chunk_0002_R2.fastq.gz
    """
}


/*
 * Gather BAM files (used when scattering is implemented)
 */
process GATHER_BAM_FILES {
    tag "${base_name}"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::picard=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/picard-cloud:2.23.8' }"
    
    input:
    path input_bams
    val base_name
    
    output:
    path "${base_name}.gathered.bam", emit: gathered_bam
    path "${base_name}.gathered.bam.bai", emit: gathered_bai
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard GatherBamFiles] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    
    def input_list = input_bams.collect { "INPUT=$it" }.join(' ')
    
    """
    java -Xmx${avail_mem}M -jar /usr/picard/picard.jar GatherBamFiles \\
        ${input_list} \\
        OUTPUT=${base_name}.gathered.bam \\
        CREATE_INDEX=false \\
        VALIDATION_STRINGENCY=SILENT \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(java -jar /usr/picard/picard.jar GatherBamFiles --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${base_name}.gathered.bam
    touch ${base_name}.gathered.bam.bai
    touch versions.yml
    """
}
