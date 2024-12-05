cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: ResourceRequirement
    ramMin: 7024
    coresMin: 1
  - class: InlineJavascriptRequirement
    expressionLib:
    - var get_target_name = function() {
          return inputs.target_filename.split('/').slice(-1)[0];
      }


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.3


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      exec 1> error_msg.txt 2>&1
      printf "rename.cwl\n$(date)\n"
      cp $0 $1
      if [ -f $0.bai ]; then
        cp $0.bai $1.bai
      fi
    inputBinding:
      position: 1

  source_file:
    type: File
    inputBinding:
      position: 5

  target_filename:
    type: string
    inputBinding:
      position: 6
      valueFrom: $(get_target_name())


outputs:

  error_msg:
    type: File?
    outputBinding:
      glob: "error_msg.txt"

  error_report:
    type: File?
    outputBinding:
      glob: "error_report.txt"

  target_file:
    type: File
    outputBinding:
      glob: $(get_target_name())
    secondaryFiles: |
      ${
          if (inputs.source_file.secondaryFiles && inputs.source_file.secondaryFiles.length > 0){
            return inputs.target_filename+".bai";
          } else {
            return "null";
          }
        }


baseCommand: [bash, '-c']


label: "rename"
doc: |
  Tool renames `source_file` to `target_filename`.
  Input `target_filename` should be set as string. If it's a full path, only basename will be used.
  If BAI file is present, it will be renamed too