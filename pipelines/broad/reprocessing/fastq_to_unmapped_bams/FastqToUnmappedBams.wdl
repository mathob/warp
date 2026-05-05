version 1.0

# Exactly one of pair of FASTQs (R1 and R2) should be supplied to this workflow.

workflow FastqToUnmappedBams {

  String pipeline_version = "1.0.0"

  input {
    File input_R1
    File input_R2
    String RGID
    String RGPL
    String RGPU
    String RGLB
    String RGCN
    String sample_name
    String zones

    String base_file_name
    String unmapped_bam_suffix = ".unmapped.bam"

    Int additional_disk = 20
  }

  String unmapped_bam_filename = base_file_name + unmapped_bam_suffix

  call FastqToSam {
    input:
      input_R1 = input_R1,
      input_R2 = input_R2,
      output_bam_filename = unmapped_bam_filename,
      RGID = RGID,
      RGLB = RGLB,
      RGPU = RGPU,
      RGPL = RGPL,
      RGCN = RGCN,
      sample_name = sample_name,
      disk_size = ceil((size(input_R1, "GiB") + size(input_R2, "GiB"))* 4) + additional_disk,
      zones = zones
  }

  Float unmapped_bam_size = size(FastqToSam.output_bam, "GiB")

  call ValidateSamFile {
    input:
      input_bam = FastqToSam.output_bam,
      report_filename = unmapped_bam_filename + ".validation_report",
      disk_size = ceil(unmapped_bam_size) + additional_disk,
      zones = zones
  }

  output {
    Array[File] validation_reports = [ValidateSamFile.report]
    Array[File] unmapped_bams = [FastqToSam.output_bam]
  }
  meta {
    allowNestedInputs: true
  }
}

task FastqToSam {
  input {
    File input_R1
    File input_R2
    String output_bam_filename
    String RGID
    String RGLB
    String RGPU
    String RGPL
    String RGCN
    String sample_name

    Int disk_size
    Int memory_in_MiB = 3000
    String zones
  }

  Int java_mem = memory_in_MiB - 1000

  command <<<
    java -Xmx6000m -jar /usr/picard/picard.jar \
      FastqToSam \
      --FASTQ ~{input_R1} \
      --FASTQ2 ~{input_R2} \
      --OUTPUT ~{output_bam_filename} \
      --READ_GROUP_NAME ~{RGID} \
      --SAMPLE_NAME ~{sample_name} \
      --LIBRARY_NAME ~{RGLB} \
      --PLATFORM_UNIT ~{RGPU} \
      --PLATFORM ~{RGPL} \
      --SEQUENCING_CENTER ~{RGCN} \
      --SORT_ORDER queryname
  >>>
  runtime {
    docker: "australia-southeast1-docker.pkg.dev/pb-dev-312200/warp/picard-cloud:2.26.10"
    disks: "local-disk " + 800 + " HDD"
    memory: "6.5 GB"
    preemptible: 3
  }

  output {
    File output_bam = output_bam_filename
  }
}

task RevertSam {
  input {
    File input_bam
    String output_bam_filename
    Int disk_size
    Int memory_in_MiB = 3000
    Int maxRetries = 1
    String zones
  }

  Int java_mem = memory_in_MiB - 1000

  command <<<
    java -Xms~{java_mem}m -jar /usr/picard/picard.jar \
    RevertSam \
    --INPUT ~{input_bam} \
    --OUTPUT ~{output_bam_filename} \
    --VALIDATION_STRINGENCY LENIENT \
    --ATTRIBUTE_TO_CLEAR FT \
    --ATTRIBUTE_TO_CLEAR CO \
    --ATTRIBUTE_TO_CLEAR PA \
    --ATTRIBUTE_TO_CLEAR OA \
    --ATTRIBUTE_TO_CLEAR XA \
    --SORT_ORDER coordinate

  >>>

  runtime {
    docker: "australia-southeast1-docker.pkg.dev/pb-dev-312200/warp/picard-cloud:2.26.10"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory_in_MiB} MiB"
    maxRetries: maxRetries
    preemptible: 3
    zones: zones
  }

  output {
    File output_bam = output_bam_filename
  }
}

# This task is slower than converting straight from cram to bam (note we stream through sam format in between cram and bam)
# This is currently necessary due to a difference in the way the NM tag is calculated in order to produce a valid bam.
task CramToBam {
  input {
    File ref_fasta
    File ref_fasta_index
    File cram_file
    String output_basename
    Int disk_size
    Int memory_in_MiB = 7000
    Int maxRetries = 1
    String zones
  }

  command <<<

    set -e
    set -o pipefail

    samtools view -h -T ~{ref_fasta} ~{cram_file} |
    samtools view -b -o ~{output_basename}.bam -
    samtools index -b ~{output_basename}.bam
    mv ~{output_basename}.bam.bai ~{output_basename}.bai

  >>>

  runtime {
    docker: "australia-southeast1-docker.pkg.dev/pb-dev-312200/warp/samtools:1.0.0-1.11-1624651616"
    cpu: 3
    memory: "~{memory_in_MiB} MiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
    maxRetries: maxRetries
    zones: zones
  }

  output {
    File output_bam = "~{output_basename}.bam"
    File output_bam_index = "~{output_basename}.bai"
  }
}

task GenerateOutputMap {
  input {
    File input_bam
    String unmapped_bam_suffix
    Int disk_size
    Int memory_in_MiB = 3000
    String zones
  }

  command <<<

    set -e

    samtools view -H ~{input_bam} | grep '^@RG' | cut -f2 | sed s/ID:// > readgroups.txt

    echo -e "#READ_GROUP_ID\tOUTPUT" > output_map.tsv

    for rg in `cat readgroups.txt`; do
      echo -e "$rg\t$rg~{unmapped_bam_suffix}" >> output_map.tsv
    done

  >>>

  runtime {
    docker: "australia-southeast1-docker.pkg.dev/pb-dev-312200/warp/samtools:1.0.0-1.11-1624651616"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory_in_MiB} MiB"
    preemptible: 3
    zones: zones
  }

  output {
    File output_map = "output_map.tsv"
  }
}

task SplitUpOutputMapFile {
  input {
    File read_group_map_file
    Int disk_size = 10
    Int memory_in_MiB = 3000
    String zones
  }

  command <<<
    mkdir rgtemp
    cd rgtemp

    # splits list of mappings into single files.  One line each.
    grep -v '^#' ~{read_group_map_file} | split -l 1 - rg_to_ubam_
  >>>

  runtime {
    docker: "australia-southeast1-docker.pkg.dev/pb-dev-312200/warp/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory_in_MiB} MiB"
    zones: zones
  }

  output {
    Array[File] rg_to_ubam_file = glob("rgtemp/rg_to_ubam_*")
  }
}

task SplitOutUbamByReadGroup {

  input {
    File input_bam
    File rg_to_ubam_file
    Int disk_size
    Int memory_in_MiB = 30000
    String zones
  }

  Array[Array[String]] tmp = read_tsv(rg_to_ubam_file)

  command <<<
    echo "Read Group ~{tmp[0][0]} from ~{input_bam} is being written to ~{tmp[0][1]}"
    samtools view -b -h -r ~{tmp[0][0]} -o ~{tmp[0][1]} ~{input_bam}
  >>>

  output {
    File output_bam = tmp[0][1]
  }

  runtime {
    docker: "australia-southeast1-docker.pkg.dev/pb-dev-312200/warp/samtools:1.0.0-1.11-1624651616"
    cpu: 2
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory_in_MiB} MiB"
    preemptible: 3
    zones: zones
  }
}

task ValidateSamFile {
  input {
    File input_bam
    String report_filename
    Int disk_size
    Int memory_in_MiB = 3000
    String zones
  }

  Int java_mem = memory_in_MiB - 1000

  command <<<

    java -Xms~{java_mem}m -jar /usr/picard/picard.jar \
      ValidateSamFile \
      --INPUT ~{input_bam} \
      --OUTPUT ~{report_filename} \
      --MODE VERBOSE \
      --IS_BISULFITE_SEQUENCED false

  >>>

  runtime {
    docker: "australia-southeast1-docker.pkg.dev/pb-dev-312200/warp/picard-cloud:2.26.10"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory_in_MiB} MiB"
    preemptible: 3
    zones: zones
    continueOnReturnCode: true
  }

  output {
    File report = "~{report_filename}"
  }
}

task SortSam {
  input {
    File input_bam
    String output_bam_filename
    Int memory_in_MiB = 7000
    Float sort_sam_disk_multiplier = 6
    Int maxRetries = 1
    String zones
  }
  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
  Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20
  Int java_mem = memory_in_MiB - 1000

  command <<<
    java -Xms~{java_mem}m -jar /usr/picard/picard.jar \
    SortSam \
    --INPUT ~{input_bam} \
    --OUTPUT ~{output_bam_filename} \
    --SORT_ORDER queryname \
    --MAX_RECORDS_IN_RAM 300000

  >>>

  runtime {
    docker: "australia-southeast1-docker.pkg.dev/pb-dev-312200/warp/picard-cloud:2.26.10"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory_in_MiB} MiB"
    preemptible: 3
    zones: zones
    maxRetries: maxRetries
  }

  output {
    File output_bam = output_bam_filename
  }
}
