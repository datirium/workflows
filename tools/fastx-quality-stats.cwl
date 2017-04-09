#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/fastx_toolkit:v0.0.14
  dockerFile: >
    $import: ./dockerfiles/fastx-Dockerfile

inputs:

  input_file:
    type: File
    inputBinding:
      position: 10
      prefix: -i
    doc: |
      FASTA/Q input file. If FASTA file is given, only nucleotides distribution is calculated (there's no quality info)

  new_output_format:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 5
      prefix: '-N'
    doc: |
      New output format (with more information per nucleotide/cycle).
      cycle (previously called 'column') = cycle number
      max-count
      For each nucleotide in the cycle (ALL/A/C/G/T/N):
          count   = number of bases found in this column.
          min     = Lowest quality score value found in this column.
          max     = Highest quality score value found in this column.
          sum     = Sum of quality score values for this column.
          mean    = Mean quality score value for this column.
          Q1	= 1st quartile quality score.
          med	= Median quality score.
          Q3	= 3rd quartile quality score.
          IQR	= Inter-Quartile range (Q3-Q1).
          lW	= 'Left-Whisker' value (for boxplotting).
          rW	= 'Right-Whisker' value (for boxplotting).

outputs:
  statistics:
    type: File
    outputBinding:
      glob: "*.fastxstat"
    doc: Statistics file

baseCommand: [fastx_quality_stats]
arguments:
  - valueFrom: $(inputs.input_file.basename + ".fastxstat")
    position: 11
    prefix: -o
