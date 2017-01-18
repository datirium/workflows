#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: envvar-global.yml

hints:
- class: DockerRequirement
  dockerPull: scidap/fastqc:v0.11.2
  dockerFile: >
    $import: fastqc-Dockerfile

inputs:
  fastqFile:
    type: File # No reason to accept multiple files as no overall report is generated
    inputBinding:
      position: 1

baseCommand: [fastqc, --outdir, ., --extract]
outputs:
  zippedFile:
    type: File
    outputBinding:
      glob: '*.zip'
  report:
    type: Directory
    outputBinding:
      glob: .