#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
    - $(inputs.input_file)

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.2
  dockerFile: >
    $import: ./dockerfiles/scidap-Dockerfile

inputs:

  input_file:
    type:
      - File
    inputBinding:
      position: 1
      valueFrom: $(self.basename)
    doc: |
      File(s) to be compressed

outputs:

  output:
    type:
      - File
    outputBinding:
      glob: |
        ${
          return inputs.input_file.basename + '.bz2';
        }

baseCommand: [bzip2, -k]