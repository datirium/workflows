#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ../../common/envvar-global.yml
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: scidap/fastx_toolkit:v0.0.14
  dockerFile: >
    $import: ../fastx-Dockerfile

inputs:

  inputFile:
    type: File
    inputBinding:
      position: 1
      prefix: -i
    doc: |
      FASTA/Q input file. If FASTA file is given, only nucleotides distribution is calculated (there's no quality info)
      If array - takes only first file, ignores all others

  outputFileDir:
    type: string
    inputBinding:
      position: 2
      prefix: -o
      valueFrom: |
        ${
          var spacer = "/";
          if (inputs.outputFileDir.slice(-1) == "/"){
            spacer = "";
          }
          return inputs.outputFileDir + spacer + inputs.inputFile.basename + ".statistics"
        }
    doc: |
      Output filename formed on the base of outputFileDir and inputFile

outputs:
  statistics:
    type: File
    outputBinding:
      glob: "*.statistics"
    doc: Statistics file

baseCommand: [fastx_quality_stats]




