cwlVersion: v1.0
class: CommandLineTool
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedops:v2.4.34
inputs:
  script:
    type: string?
    default: |
      cat "$0" > `basename $0`
    inputBinding:
      position: 1
  input_file:
    type:
    - File
    - File[]
    inputBinding:
      position: 2
  param:
    type:
    - string?
    - string[]
    inputBinding:
      position: 3
outputs:
  output_file:
    type: File
    outputBinding:
      glob: '*'
baseCommand:
- bash
- -c
doc: |
  Tool to run custom script set as `script`
  input with arguments from `param`. Based
  on bedops Dockerfile.
label: custom-bedops
