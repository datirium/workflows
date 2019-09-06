cwlVersion: v1.0
class: ExpressionTool


requirements:
  - class: InlineJavascriptRequirement


inputs:

  input_files:
    type:
      - Directory
      - File[]

  basename_regex:
    type: string


outputs:

  selected_file:
    type: File?


expression: |
  ${
    var patt = new RegExp(inputs.basename_regex);
    var files = [];

    if (inputs.input_files.class == "Directory"){
      files = inputs.input_files.listing;
    } else {
      files = inputs.input_files;
    }

    for (var i = 0; i < files.length; i++ ){
      if ( patt.test(files[i].location.split('/').slice(-1)[0]) ){
        return { "selected_file": files[i] }
      }
    }

    return { "selected_file": null };
  }


doc: |
  Returns file the first file from input File[], that match input regex expression