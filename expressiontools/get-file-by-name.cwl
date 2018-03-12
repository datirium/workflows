cwlVersion: v1.0
class: ExpressionTool

requirements:
  - class: InlineJavascriptRequirement

class: ExpressionTool
id: "get_file_by_name"
inputs:
  input_files: File[]
  basename_regex: string
outputs:
  selected_file: File
expression: |
  ${
    var patt = new RegExp(inputs.basename_regex);
    for (var i = 0; i < inputs.input_files.length; i++ ){
      if ( patt.test(inputs.input_files[i].location.split('/').slice(-1)[0]) ){
        return { "selected_file": inputs.input_files[i] }
      }
    }
    return null
  }

doc: |
  Returns file the first file from input File[], that match input regex expression