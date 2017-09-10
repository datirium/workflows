cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement


inputs:

  sra_input_file:
    type: File
  illumina_adapters_file:
    type: File
  bowtie2_indices_folder:
    type: Directory
  chr_length_file:
    type: File
  threads:
    type: int?

outputs:

  fastq:
    type: File
    outputSource: sra_to_fastq/fastq
  bowtie2_log:
    type: File
    outputSource: fastq_to_bigwig/bowtie2_log
  picard_metrics:
    type: File
    outputSource: fastq_to_bigwig/picard_metrics
  bam_file:
    type: File
    outputSource: fastq_to_bigwig/bam_file
  bamtools_stats_log:
    type: File
    outputSource: fastq_to_bigwig/bamtools_stats_log
  bed:
    type: File
    outputSource: fastq_to_bigwig/bed
  bigwig:
    type: File
    outputSource: fastq_to_bigwig/bigwig


steps:

  sra_to_fastq:
    run: xenbase-sra-to-fastq-se.cwl
    in:
      sra_input_file: sra_input_file
      illumina_adapters_file: illumina_adapters_file
      threads: threads
    out: [fastq]

  fastq_to_bigwig:
    run: xenbase-fastq-bowtie-bigwig-se-pe.cwl
    in:
      upstream_fastq: sra_to_fastq/fastq
      bowtie2_indices_folder: bowtie2_indices_folder
      chr_length_file: chr_length_file
      threads: threads
    out:
      - bowtie2_log
      - picard_metrics
      - bam_file
      - bamtools_stats_log
      - bed
      - bigwig