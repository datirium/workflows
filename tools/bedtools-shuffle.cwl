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
      The input BED file with features to be shuffled
  chrom_length_file:
    type: File
    inputBinding:
      position: 6
      prefix: -g
    doc: |
      Input genome file with chromosome lengths
  excl_bed_file:
    type: File?
    inputBinding:
      position: 7
      prefix: -excl
    doc: |
      The BED file with regions where you do not want the permuted features to be placed
  incl_bed_file:
    type: File?
    inputBinding:
      position: 8
      prefix: -incl
    doc: |
      The BED file with regions where you want the permuted features to be placed
  seed:
    type: int?
    inputBinding:
      position: 9
      prefix: -seed
    doc: |
      Seed for pseudo-random number generation. Default: random
  no_overlapping:
    type: boolean?
    inputBinding:
      position: 10
      prefix: -noOverlapping
    doc: |
      Don't allow shuffled intervals to overlap
  max_tries:
    type: int?
    inputBinding:
      position: 11
      prefix: -maxTries
    doc: |
      Maximum number of attempts to find a home for a shuffled interval. Default: 1000
  output_filename:
    type: string?
    default: ''
    doc: |
      Output file name
outputs:
  shuffled_bed_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: |
      Shuffled BED file
baseCommand:
- bedtools
- shuffle
stdout: $(default_output_filename())
doc: |
  Randomly permutes the genomic locations of a feature file among a genome defined in a genome file.
  One can also provide an “exclusions” BED file that lists regions where you do not want the permuted
  features to be placed. Or instead “inclusions” BED fils that defines coordinates in which features
  in -i should be randomly placed. To make experiment reproducible, set "seed" option.
  NOTE: limited parameters are impelented.
label: bedtools-shuffle
