cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() { var basename = inputs.bigwig_file.location.split('/').slice(-1)[0]; var root = basename.split('.').slice(0,-1).join('.'); var ext = ".wig"; return (root == "")?basename+ext:root+ext; };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/ucscuserapps:v358
inputs:
  bigwig_file:
    type: File
    inputBinding:
      position: 1
    doc: |
      Input bigWig file
  chrom:
    type: string?
    inputBinding:
      position: 2
      prefix: -chrom=
      separate: false
    doc: |
      if set restrict output to given chromosome
  start_pos:
    type: int?
    inputBinding:
      position: 3
      prefix: -start=
      separate: false
    doc: |
      if set, restrict output to only that over start
  end_pos:
    type: int?
    inputBinding:
      position: 4
      prefix: -end=
      separate: false
    doc: |
      if set, restict output to only that under end
  output_filename:
    type: string?
    inputBinding:
      position: 5
      valueFrom: |
        ${
            if (self == ""){
              return default_output_filename();
            } else {
              return self;
            }
        }
    default: ''
    doc: |
      If set, writes the output wig file to output_filename,
      otherwise generates filename from default_output_filename()
outputs:
  wig_file:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == ""){
              return default_output_filename();
            } else {
              return inputs.output_filename;
            }
        }
baseCommand:
- bigWigToWig
doc: |
  Tool converts bigWig to Wig file. If input bigWig file was generated from bedGraph, the tool will
  return output in bedGraph format.

  `default_output_filename` function returns filename for generated Wig if `output_filename` is not provided.
  Default filename is generated on the base of `bigwig_file` basename with the updated to `*.wig` extension.
label: ucsc-bigwigtowig
