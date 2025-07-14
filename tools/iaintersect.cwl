cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function(ext) { var root = inputs.input_filename.basename.split('.').slice(0,-1).join('.'); return (root == "")?inputs.input_filename.basename+ext:root+ext; };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/iaintersect:v0.0.2
inputs:
  input_filename:
    type: File
    inputBinding:
      position: 1
      prefix: --in=
      separate: false
    doc: Input filename with MACS2 peak calling results, tsv
  annotation_filename:
    type: File
    inputBinding:
      position: 2
      prefix: --a=
      separate: false
    doc: Annotation file, tsv
  output_filename:
    type: string?
    inputBinding:
      position: 3
      prefix: --out=
      separate: false
      valueFrom: |
        ${
            if (self == ""){
              return default_output_filename('_iaintersect.tsv');
            } else {
              return self;
            }
        }
    default: ''
    doc: |
      Base output file name, tsv
  log_filename:
    type: string?
    inputBinding:
      position: 4
      prefix: --log=
      separate: false
      valueFrom: |
        ${
            if (self == ""){
              return default_output_filename('_iaintersect.log');
            } else {
              return self;
            }
        }
    default: ''
    doc: |
      Log filename
  promoter_bp:
    type: int?
    inputBinding:
      position: 5
      prefix: --promoter=
      separate: false
    doc: |
      Promoter region around TSS, base pairs
  upstream_bp:
    type: int?
    inputBinding:
      position: 6
      prefix: --upstream=
      separate: false
    doc: |
      Upstream region before promoter, base pairs
  ignore_chr:
    type: string?
    inputBinding:
      position: 7
      prefix: --sam_ignorechr=
      separate: false
    doc: |
      The chromosome to be ignored, string
outputs:
  log_file:
    type: File?
    outputBinding:
      glob: |
        ${
          if (inputs.log_filename == ""){
            return default_output_filename('_iaintersect.log');
          } else {
            return inputs.log_filename;
          }
        }
  result_file:
    type: File?
    outputBinding:
      glob: |
        ${
          if (inputs.output_filename == ""){
            return default_output_filename('_iaintersect.tsv');
          } else {
            return inputs.output_filename;
          }
        }
baseCommand:
- iaintersect
doc: |
  Tool assigns each peak obtained from MACS2 to a gene and region (upstream, promoter, exon, intron, intergenic)

  `default_output_filename` function returns output filename with sufix set as `ext` argument. Function is called when
  either `output_filename` or `log_filename` inputs are not provided.
label: iaintersect
