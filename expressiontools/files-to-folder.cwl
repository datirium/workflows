cwlVersion: v1.0
class: ExpressionTool

requirements:
  - class: InlineJavascriptRequirement

class: ExpressionTool
id: "files_to_folder"
inputs:
  input_files: File[]
outputs:
  folder: Directory
expression: |
  ${
    var folder = {
      "class": "Directory",
      "basename": "folder",
      "listing": []
    }
    folder.listing = inputs.input_files;
    return { "folder": folder };
  }