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
  bowtie2_indices:
    type: File
    doc: "Path to any of the BOWTIE2 index file"
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
  output_peak_file_narrow:
    type: File
    outputSource: macs2_callpeak_narrow/output_peak_file
  output_peak_file_broad:
    type: File
    outputSource: macs2_callpeak_broad/output_peak_file


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
      bowtie2_indices: bowtie2_indices
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

  macs2_callpeak_narrow:
    run: ../../tools/macs2-callpeak.cwl
    in:
      treatment: from_bowtie2_to_bigwig_SE_PE/bam_file
      format:
        default: BAM
    out: [output_peak_file, output_peak_xls_file, output_peak_summits_file]

  macs2_callpeak_broad:
    run: ../../tools/macs2-callpeak.cwl
    in:
      treatment: from_bowtie2_to_bigwig_SE_PE/bam_file
      broad:
        default: true
      format:
        default: BAM
    out: [output_peak_file, output_peak_xls_file, output_peak_summits_file]