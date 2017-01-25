#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: envvar-global.yml
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: scidap/parse-fastqc-results:v0.0.1
  dockerFile: >
    $import: parse-fastqc-results-Dockerfile

inputs:
  fastqc_report:
    type: File
    inputBinding:
      position: 1
    doc: >
      Path to the summary.txt file, where fastqc exported its results

baseCommand: [python, /usr/local/bin/fastqc_pass_test_trigger.py]

outputs:
  output:
    type: File
    outputBinding:
      glob: "*.txt"