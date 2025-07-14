cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ResourceRequirement
  ramMin: 7620
  coresMin: 1
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_filename = function() { if (inputs.output_filename) { return inputs.output_filename; } return inputs.unsorted_file.location.split('/').slice(-1)[0]; };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.2
inputs:
  unsorted_file:
    type: File
    inputBinding:
      position: 4
  key:
    type:
      type: array
      items: string
      inputBinding:
        prefix: -k
    inputBinding:
      position: 1
    doc: |
      -k, --key=POS1[,POS2]
      start a key at POS1, end it at POS2 (origin 1)
  output_filename:
    type: string?
    doc: |
      Name for generated output file
outputs:
  sorted_file:
    type: stdout
stdout: $(get_output_filename())
baseCommand:
- sort
doc: |
  Tool sorts data from `unsorted_file` by key

  `default_output_filename` function returns file name identical to `unsorted_file`, if `output_filename` is not provided.
label: linux-sort
