cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.3
- class: InlineJavascriptRequirement
- class: ResourceRequirement
  ramMin: $(inputs.memory)
  coresMin: $(inputs.threads)

inputs:

  delay:
    type: int?
    default: 60
    inputBinding:
      position: 1
    doc: "Delay"

  threads:
    type: int?
    default: 1
    doc: "CPU"

  memory:
    type: int?
    default: 15259             # equal to 16GB
    doc: "Memory"


outputs:

  dummy:
    type: File?
    outputBinding:
      glob: "*"

baseCommand: [sleep]
