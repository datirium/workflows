#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: envvar-global.yml
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: scidap/fastqc:v0.11.2
  dockerFile: >
    $import: ./dockerfiles/fastqc-Dockerfile

inputs:
  fastqFile:
    type: File
    inputBinding:
      position: 1

baseCommand: [fastqc, --outdir, ., --extract]
outputs:
  zippedFile:
    type: File
    outputBinding:
      glob: '*.zip'
#  report_dir:
#    type: Directory
#    outputBinding:
#      glob: .
  summary_file:
    type: File
    outputBinding:
      glob: |
        ${
          return inputs.fastqFile.nameroot + "_fastqc/summary.txt";
        }