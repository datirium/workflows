cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() { if (inputs.output_filename == ""){ return inputs.file_a.basename; } else { return inputs.output_filename; } };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0
inputs:
  file_a:
    type: File
    inputBinding:
      position: 5
      prefix: -a
    doc: BAM/BED/GFF/VCF file A. Each feature in A is compared to B in search of overlaps
  file_b:
    type: File
    inputBinding:
      position: 6
      prefix: -b
    doc: BAM/BED/GFF/VCF file B. Each feature in A is compared to B in search of overlaps
  count:
    type: boolean?
    inputBinding:
      position: 7
      prefix: -c
    doc: For each entry in A, report the number of hits in B. Reports 0 for A entries that have no overlap with B
  no_overlaps:
    type: boolean?
    inputBinding:
      position: 8
      prefix: -v
    doc: Only report those entries in A that have _no overlaps_ with B
  report_from_a_once:
    type: boolean?
    inputBinding:
      position: 9
      prefix: -u
    doc: Write the original A entry _once_ if _any_ overlaps found in B
  output_filename:
    type: string?
    default: ''
    doc: Output file name
outputs:
  intersected_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: Intersected BED file
baseCommand:
- bedtools
- intersect
stdout: $(default_output_filename())
doc: |
  Intersect features from A and B file. Only selected parameters are implemented.
label: bedtools-intersect
