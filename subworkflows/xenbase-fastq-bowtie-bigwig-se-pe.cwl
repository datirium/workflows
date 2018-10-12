cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement


inputs:
  upstream_fastq:
    type: File
  downstream_fastq:
    type: File?
  bowtie2_indices_folder:
    type: Directory
  chr_length_file:
    type: File
  paired:
    type: boolean?
  threads:
    type: int?


outputs:

  bowtie2_log:
    type: File
    outputSource: bowtie2_aligner/output_log
  picard_metrics:
    type: File
    outputSource: remove_dup_picard/metrics_file
  bam_file:
    type: File
    outputSource: samtools_sort_index_after_dup_removing/bam_bai_pair
  bamtools_log:
    type: File
    outputSource: bamtools_stats/log_file
  bed:
    type: File
    outputSource: bam_to_bigwig/bedgraph_file
  bigwig:
    type: File
    outputSource: bam_to_bigwig/bigwig_file


steps:
  bowtie2_aligner:
    run: ../tools/bowtie2.cwl
    in:
      filelist: upstream_fastq
      filelist_mates: downstream_fastq
      indices_folder: bowtie2_indices_folder
      threads: threads
    out: [output, output_log]

  samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: bowtie2_aligner/output
      threads: threads
    out: [bam_bai_pair]

  remove_dup_picard:
    run: ../tools/picard-markduplicates.cwl
    in:
      input_file: samtools_sort_index/bam_bai_pair
      remove_dup:
        default: true
    out:
      - output_file
      - metrics_file

  samtools_sort_index_after_dup_removing:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: remove_dup_picard/output_file
      threads: threads
    out: [bam_bai_pair]

  bamtools_stats:
    run: ../tools/bamtools-stats.cwl
    in:
      bam_file: samtools_sort_index_after_dup_removing/bam_bai_pair
    out:
      - mapped_reads_number
      - log_file

  bam_to_bigwig:
    run: bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index_after_dup_removing/bam_bai_pair
      chrom_length_file: chr_length_file
      mapped_reads_number: bamtools_stats/mapped_reads_number
      pairchip: paired
    out:
      - bigwig_file
      - bedgraph_file