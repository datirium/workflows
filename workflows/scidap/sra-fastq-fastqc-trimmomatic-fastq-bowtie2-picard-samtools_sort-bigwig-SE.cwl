cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  sra_input_file:
    type: File
  illumina_adapters_file:
    type: File
  bowtie2_indices_folder:
    type: Directory
    doc: "Path to BOWTIE2 indices folder"
  samtools_view_reads_quality_cutoff:
    type: int
    doc: Filter out all aligned reads, whch have quality lower then this threshold
  chrLengthFile:
    type: File

outputs:
  fastq:
    type: File
    outputSource: sra_fastqc_trimmomatic_fastq_SE/fastq
  bowtie2_aligner_log:
    type: File
    outputSource: from_bowtie2_to_bigwig_SE_PE/bowtie2_aligner_log
  picard_metrics:
    type: File
    outputSource: from_bowtie2_to_bigwig_SE_PE/picard_metrics
  bam_file:
    type: File
    outputSource: from_bowtie2_to_bigwig_SE_PE/bam_file
  bamtools_stats_log:
    type: File
    outputSource: from_bowtie2_to_bigwig_SE_PE/bamtools_stats_log
  bed:
    type: File
    outputSource: from_bowtie2_to_bigwig_SE_PE/bed
  bigwig:
    type: File
    outputSource: from_bowtie2_to_bigwig_SE_PE/bigwig

steps:

  sra_fastqc_trimmomatic_fastq_SE:
    run: sra-fastq-fastqc-trimmomatic-fastq-SE.cwl
    in:
      sra_input_file: sra_input_file
      illumina_adapters_file: illumina_adapters_file
    out: [fastq]

  from_bowtie2_to_bigwig_SE_PE:
    run: fastq-bowtie2-picard-samtools_sort-bigwig-SE-PE.cwl
    in:
      upstream_fastq: sra_fastqc_trimmomatic_fastq_SE/fastq
      bowtie2_indices_folder: bowtie2_indices_folder
      samtools_view_reads_quality_cutoff: samtools_view_reads_quality_cutoff
      chrLengthFile: chrLengthFile
      split:
        default: true
    out:
    - bowtie2_aligner_log
    - picard_metrics
    - bam_file
    - bamtools_stats_log
    - bed
    - bigwig