cwlVersion: v1.0
class: Workflow


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
    outputSource: bam_to_bigwig/outfile


  bam_to_bigwig:
    run: bam-genomecov-bigwig.cwl
    in:
      input: bam_bai_pair_file
      genomeFile: chrom_length_file
      mappedreads: mapped_reads
      fragmentsize: fragment_size
      pairchip: paired
    out: [outfile]