cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: ubuntu:20.04
inputs:
  folder_to_compress:
    type: Directory
    doc: Folder to compressed
outputs:
  compressed_folder:
    type: File
    outputBinding:
      glob: '*'
    doc: Compressed folder
baseCommand:
- tar
arguments:
- valueFrom: $(inputs.folder_to_compress.path.split("/").slice(0,-1).join("/"))
  prefix: -C
- -czvf
- valueFrom: $(inputs.folder_to_compress.basename + ".tar.gz")
- .
label: TAR compress
doc: |
  TAR compress
  =========================================

  Creates compressed TAR file from a folder
