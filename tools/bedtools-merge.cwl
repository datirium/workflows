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
    doc: The input BED file must be sorted by chrom, then start
  max_distance:
    type: int?
    inputBinding:
      position: 6
      prefix: -d
    doc: Maximum distance between features to be merged
  output_filename:
    type: string?
    default: ''
    doc: Output file name
outputs:
  merged_bed_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: Merged BED file
baseCommand:
- bedtools
- merge
stdout: $(default_output_filename())
doc: |
  Merges features from BED file. Only selected parameters are implemented.
label: bedtools-merge
