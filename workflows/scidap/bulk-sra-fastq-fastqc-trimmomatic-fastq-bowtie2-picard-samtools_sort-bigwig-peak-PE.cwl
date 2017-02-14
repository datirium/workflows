cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  sra_input_files:
    type: File[]
  illumina_adapters_file:
    type: File
  bowtie2_indices:
    type: File
    doc: "Path to any of the BOWTIE2 index file"
  samtools_view_reads_quality_cutoff:
    type: int
    doc: Filter out all aligned reads, whch have quality lower then this threshold
  chrLengthFile:
    type: File

outputs:
  upstream_fastq:
    type: File[]
    outputSource: subworkflow-PE/upstream_fastq
  downstream_fastq:
    type: File[]
    outputSource: subworkflow-PE/downstream_fastq
  bowtie2_aligner_log:
    type: File[]
    outputSource: subworkflow-PE/bowtie2_aligner_log
  picard_metrics:
    type: File[]
    outputSource: subworkflow-PE/picard_metrics
  bam_file:
    type: File[]
    outputSource: subworkflow-PE/bam_file
  bamtools_stats_log:
    type: File[]
    outputSource: subworkflow-PE/bamtools_stats_log
  bed:
    type: File[]
    outputSource: subworkflow-PE/bed
  bigwig:
    type: File[]
    outputSource: subworkflow-PE/bigwig
  output_peak_file_narrow:
    type: File[]
    outputSource: subworkflow-PE/output_peak_file_narrow
  output_peak_file_broad:
    type: File[]
    outputSource: subworkflow-PE/output_peak_file_broad


steps:

  subworkflow-PE:
    run: sra-fastq-fastqc-trimmomatic-fastq-bowtie2-picard-samtools_sort-bigwig-peak-PE.cwl
    in:
      sra_input_file: sra_input_files
      illumina_adapters_file: illumina_adapters_file
      bowtie2_indices: bowtie2_indices
      samtools_view_reads_quality_cutoff: samtools_view_reads_quality_cutoff
      chrLengthFile: chrLengthFile
    scatter: sra_input_file
    out:
    - upstream_fastq
    - downstream_fastq
    - bowtie2_aligner_log
    - picard_metrics
    - bam_file
    - bamtools_stats_log
    - bed
    - bigwig
    - output_peak_file_narrow
    - output_peak_file_broad