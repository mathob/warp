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
    samtools view -@ ${task.cpus} -Sb - > ${sample_name}.aligned.bam
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_name}.aligned.bam
    touch versions.yml
    """
}

/*
 * DRAGMAP alignment (alternative to BWA MEM)
 */
process DRAGMAP_ALIGN {
    tag "${sample_name}"
    #label 'process_high'
    
    conda (params.enable_conda ? "bioconda::dragmap=1.2.1 bioconda::samtools=1.17" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/dragmap:1.2.1' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/dragmap:1.2.1' }"
    
    input:
    path unmapped_bam
    path dragmap_reference_bin
    path dragmap_hash_table_cfg_bin
    path dragmap_hash_table_cmp
    path reference_fasta
    path reference_dict
    val sample_name
    val read_group_id
    val read_group_platform
    val read_group_pu
    val read_group_library
    val read_group_center
    
    output:
    path "${sample_name}.chunk_*.aligned.bam", emit: aligned_bam
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def chunk_id   = (unmapped_bam.name =~ /chunk_(\d+).unmapped/)[0]?[1] ?: "000"
    def aligned_unmerged_bam = "${sample_name}.chunk_${chunk_id}.aligned.unmerged.bam"
    def output_bam = "${sample_name}.chunk_${chunk_id}.aligned.bam"
    //def read_group = "@RG\\tID:${read_group_id}\\tSM:${sample_name}\\tPL:${read_group_platform}\\tPU:${read_group_pu}\\tLB:${read_group_library}\\tCN:${read_group_center}"
    avail_mem = (task.memory.mega*0.8).intValue()

    """
    # DRAGMAP alignment
    dragen-os \\
        -b ${unmapped_bam} \\
        --ref-dir . \\
        --RGID ${read_group_id} \\
        --RGSM ${sample_name} \\
        --interleaved 1 \\
        --preserve-map-align-order true \\
        --num-threads ${task.cpus} \\
        ${args} \\
    2> ${sample_name}.dragmap.log \\
        | samtools view --threads ${task.cpus} -o ${aligned_unmerged_bam} -
    
    # Merge unmapped and aligned bams
    java -Dsamjdk.compression_level=2 -Xmx${avail_mem}M -Xms${avail_mem}M -jar /picard/picard.jar \
      MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ATTRIBUTES_TO_REMOVE=RG \
      ATTRIBUTES_TO_REMOVE=NM \
      ATTRIBUTES_TO_REMOVE=MD \
      ALIGNED_BAM=${aligned_unmerged_bam} \
      UNMAPPED_BAM=${unmapped_bam} \
      OUTPUT=${output_bam} \
      REFERENCE_SEQUENCE=${reference_fasta} \
      PAIRED_RUN=true \
      SORT_ORDER="unsorted" \
      IS_BISULFITE_SEQUENCE=false \
      ALIGNED_READS_ONLY=false \
      CLIP_ADAPTERS=false \
      MAX_RECORDS_IN_RAM=2000000 \
      ADD_MATE_CIGAR=true \
      MAX_INSERTIONS_OR_DELETIONS=-1 \
      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
      PROGRAM_RECORD_ID="dragen-os" \
      PROGRAM_GROUP_VERSION="UNKNOWN" \
      PROGRAM_GROUP_COMMAND_LINE="dragen-os -b ${unmapped_bam} -r dragen_reference --interleaved=1" \
      PROGRAM_GROUP_NAME="dragen-os" \
      UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
      ALIGNER_PROPER_PAIR_FLAGS=true \
      UNMAP_CONTAMINANT_READS=false \
      ADD_PG_TAG_TO_READS=false




    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragmap: UNKNOWN
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        picard: \$(java -jar /usr/picard/picard.jar MergeBamAlignment --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_name}.aligned.bam
    touch versions.yml
    """
}


/*
 * Calibrate DragSTR model
 * Equivalent to CalibrateDragstrModel in tasks/broad/DragenTasks.wdl
 */
process CALIBRATE_DRAGSTR_MODEL {

    tag "${sample_name}"

    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.4.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.2.2.0' :
        'australia-southeast1-docker.pkg.dev/pb-dev-312200/nagim-images/gatk:4.2.2.0' }"
    
    input:
    path bam
    path bam_index
    path reference_fasta
    path reference_fasta_index
    path reference_dict
    path str_table_file
    val  sample_name

    output:
    path "${sample_name}.dragstr_model.txt", emit: dragstr_model

    script:
    def mem_gb = task.memory ? (task.memory.toGiga() - 1) : 3
    """
    gatk --java-options "-Xmx${mem_gb}g" \
        CalibrateDragstrModel \
        -R ${reference_fasta} \
        -I ${bam} \
        -str ${str_table_file} \
        -O ${sample_name}.dragstr_model.txt \
        --parallel

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_name}.dragstr_model.txt
    """
}