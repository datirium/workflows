cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() { if (inputs.output_filename == ""){ return inputs.reduced_bed_file.basename; } else { return inputs.output_filename; } };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0
inputs:
  reduced_bed_file:
    type: File
    inputBinding:
      position: 5
      prefix: -a
    doc: |
      The input BED file from which the features will be substructed
  subtracted_bed_file:
    type: File
    inputBinding:
      position: 6
      prefix: -b
    doc: |
      The input BED file that includes the features to be subtracted from -a
  output_filename:
    type: string?
    default: ''
    doc: |
      Output file name
outputs:
  difference_bed_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: |
      Difference BED file
baseCommand:
- bedtools
- subtract
stdout: $(default_output_filename())
doc: |
  Searches for features in B that overlap A by at least 1 base pair.
  If an overlapping feature is found in B, the overlapping portion is removed from A
  and the remaining portion of A is reported. If a feature in B overlaps all of
  a feature in A, the A feature will not be reported.
  All parameters except -a and -b used by default.
label: bedtools-subtract
