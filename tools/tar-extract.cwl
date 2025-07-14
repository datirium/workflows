cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: ubuntu:20.04
inputs:
  file_to_extract:
    type: File
    inputBinding:
      position: 1
    doc: File to extract
outputs:
  extracted_folder:
    type: Directory
    outputBinding:
      glob: '*'
    doc: Extracted folder
baseCommand:
- tar
- xzf
label: TAR extract
doc: |-
  TAR extract
  ===============================================

  Extracts the content of TAR file into a folder.
