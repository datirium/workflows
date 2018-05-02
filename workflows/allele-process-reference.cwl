cwlVersion: v1.0
class: Workflow


requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement
- class: InlineJavascriptRequirement


inputs:

  bam_file:
    type: File
    label: "BAM file"
    doc: "BAM file mapped to reference genome"

  chrom_length_file:
    type: File
    label: "Chromosome length file for reference genome"
    doc: "Chromosome length file for reference genome"

  mapped_reads_number:
    type: int
    label: "Uniquely mapped reads number"
    doc: "Uniquely mapped to reference genome reads number"

  output_file_prefix:
    type: string
    label: "Prefix for all generated output files"
    doc: "Corresponds to UID"

  threads:
    type: int?
    default: 2
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"


outputs:

  bambai_pair:
    type: File
    outputSource: samtools_sort_index/bam_bai_pair
    label: "Reference BAM"
    doc: "Coordinate sorted BAM file mapped to reference genome"

  bigwig_file:
    type: File
    label: "Reference bigWig file"
    doc: "Generated bigWig file for reference genome"
    outputSource: bam_to_bigwig/bigwig_file


steps:

  samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: bam_file
      sort_output_filename:
        source: output_file_prefix
        valueFrom: $(self+".bam")
      threads: threads
    out: [bam_bai_pair]

  bam_to_bigwig:
    run: bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index/bam_bai_pair
      chrom_length_file: chrom_length_file
      mapped_reads_number: mapped_reads_number
    out: [bigwig_file]