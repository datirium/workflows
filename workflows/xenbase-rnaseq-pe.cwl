cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement


inputs:
  sra_input_file:
    type: File
  illumina_adapters_file:
    type: File
  rsem_indices_folder:
    type: Directory
  chr_length_file:
    type: File
  threads:
    type: int?


outputs:
  upstream_fastq:
    type: File
    outputSource: sra_to_fastq/upstream_fastq
  downstream_fastq:
    type: File
    outputSource: sra_to_fastq/downstream_fastq
  rsem_isoforms:
    type: File
    outputSource: fastq_to_bigwig/rsem_isoforms
  rsem_genes:
    type: File
    outputSource: fastq_to_bigwig/rsem_genes
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
    run: xenbase-sra-to-fastq-pe.cwl
    in:
      sra_input_file: sra_input_file
      illumina_adapters_file: illumina_adapters_file
      threads: threads
    out:
      - upstream_fastq
      - downstream_fastq

  fastq_to_bigwig:
    run: xenbase-fastq-rsem-bigwig-se-pe.cwl
    in:
      upstream_fastq: sra_to_fastq/upstream_fastq
      downstream_fastq: sra_to_fastq/downstream_fastq
      rsem_indices_folder: rsem_indices_folder
      chr_length_file: chr_length_file
      threads: threads
    out:
      - rsem_isoforms
      - rsem_genes
      - bam_file
      - bamtools_stats_log
      - bed
      - bigwig
