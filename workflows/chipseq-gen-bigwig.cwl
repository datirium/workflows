cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement

inputs:
  bam_bai_pair_file:
    type: File
    doc: "BAM file with index BAI file"

  chrom_length_file:
    type: File
    doc: "Chromosome length file"

  mapped_reads:
    type: int
    doc: "Mapped reads number, tagsmapped"

  fragment_size:
    type: int
    doc: "Fragment size, "

  paired:
    type: boolean
    doc: "True if paired end experiment, fragmentsize"


outputs:
  bigwig_file:
    type: File
    doc: "bigWig generated file"
    outputSource: bam_to_bigwig/bigwig_file

steps:
  bam_to_bigwig:
    run: bam-bedgraph-bigwig.cwl
    in:
      bam_file: bam_bai_pair_file
      chrom_length_file: chrom_length_file
      mapped_reads_number: mapped_reads
      fragment_size: fragment_size
      pairchip: paired
    out: [bigwig_file]