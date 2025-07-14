cwlVersion: v1.0
class: CommandLineTool
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.3
inputs:
  script:
    type: string?
    default: |
      cat "$0" | grep "$1" | sed "s/$1//g"  > `basename $0`
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
  Tool to run custom script set as `script` input with arguments from `param`.
  Default script runs sed command over the input file and exports results to the file with the same name as input's basename
label: custom-bash
