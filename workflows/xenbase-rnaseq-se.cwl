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

  rsem_isoforms:
    type: File
    outputSource: fastq_to_bigwig/rsem_isoforms
  rsem_genes:
    type: File
    outputSource: fastq_to_bigwig/rsem_genes
  bam_file:
    type: File
    outputSource: fastq_to_bigwig/bam_file
  bamtools_log:
    type: File
    outputSource: fastq_to_bigwig/bamtools_log
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
    run: xenbase-fastq-rsem-bigwig-se-pe.cwl
    in:
      upstream_fastq: sra_to_fastq/fastq
      rsem_indices_folder: rsem_indices_folder
      chr_length_file: chr_length_file
      threads: threads
    out:
      - rsem_isoforms
      - rsem_genes
      - bam_file
      - bamtools_log
      - bed
      - bigwig
