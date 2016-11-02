cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:

  genomeFastaFiles:
    type: File[]
    format:
      - http://edamontology.org/format_1929 # FASTA
    inputBinding:
      position: 1
      itemSeparator: ' '
      prefix: --genomeFastaFiles
    doc: |
      string(s): path(s) to the fasta files with genomic sequences for genome
      generation, separated by spaces. Only used if runMode==genomeGenerate.
      These files should be plain text FASTA files, they *cannot* be zipped.

  sjdbGTFfile:
    type: File?
    format:
      - http://edamontology.org/format_2306
    inputBinding:
      position: 1
      prefix: --sjdbGTFfile
    doc: |
      string: path to the GTF file with annotations

outputs:
  outfileBigWig:
    type: File[]
    outputSource: indexGenome/indices

steps:
  indexGenome:
    run: ../../tools/star-genomegenerate.cwl
    in:
      genomeFastaFiles: genomeFastaFiles
      sjdbGTFfile: sjdbGTFfile
    out: [indices]