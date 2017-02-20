cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  fastqFiles:
    type: File[]

outputs:
  zippedFiles:
    type: File[]
    outputSource: fastqc/zippedFile
  summary_files:
    type: File[]
    outputSource: fastqc/summary_file

steps:

  fastqc:
    run: ../../tools/fastqc.cwl
    in:
      fastqFile: fastqFiles
    scatter: fastqFile
    out:
    - zippedFile
    - summary_file

