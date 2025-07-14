cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() { if (inputs.output_filename == ""){ return inputs.bed_file.basename; } else { return inputs.output_filename; } };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0
inputs:
  bed_file:
    type: File
    inputBinding:
      position: 5
      prefix: -i
    doc: |
      The input BED file
  chrom_length_file:
    type: File
    inputBinding:
      position: 6
      prefix: -g
    doc: |
      Input genome file with chromosome lengths
  bi_direction:
    type:
    - 'null'
    - int
    - float
    inputBinding:
      position: 7
      prefix: -b
    doc: |
      Increase the BED entry by the same number base pairs
      or a fraction of the feature's length in each direction
  left_direction:
    type:
    - 'null'
    - int
    - float
    inputBinding:
      position: 8
      prefix: -l
    doc: |
      The number of base pairs or a fraction of the feature's length
      to subtract from the start coordinate
  right_direction:
    type:
    - 'null'
    - int
    - float
    inputBinding:
      position: 9
      prefix: -r
    doc: |
      The number of base pairs or a fraction of the feature's length
      to add to the end coordinate
  strand_based:
    type: boolean?
    inputBinding:
      position: 10
      prefix: -s
    doc: "Define -l and -r based on strand\nE.g. if used, -l 500 for a negative-stranded feature, \nit will add 500 bp downstream\n"
  percent_based:
    type: boolean?
    inputBinding:
      position: 11
      prefix: -pct
    doc: |
      Define -l, -b and -r as a fraction of the feature's length
  save_header:
    type: boolean?
    inputBinding:
      position: 12
      prefix: -header
    doc: |
      Print the header from the input file prior to results
  output_filename:
    type: string?
    default: ''
    doc: Output file name
outputs:
  extended_bed_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: Extended BED file
baseCommand:
- bedtools
- slop
stdout: $(default_output_filename())
doc: |
  Increases the size of each feature in a feature file by a user-defined number of bases.
  If not using -b, then -l and -r should be used together
label: bedtools-slop
