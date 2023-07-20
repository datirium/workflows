cwlVersion: v1.0
class: ExpressionTool

requirements:
  - class: InlineJavascriptRequirement

id: "files_to_folder"
inputs:
  input_files:
    type:
    - File[]
    - File
  folder_basename:
    type: string?
    default: ""
outputs:
  folder: Directory
expression: |
  ${
    var folder_basename = inputs.folder_basename.split('/').slice(-1).join('');
    var folder = {
      "class": "Directory",
      "basename": folder_basename,
      "listing": []
    }
    var files = [];
    if (!Array.isArray(inputs.input_files)){
      files.push (inputs.input_files)
    } else {
      files = inputs.input_files
    }
    for (var i = 0; i < files.length; i++){
      if (files[i].secondaryFiles && files[i].secondaryFiles.length > 0){
        Array.prototype.push.apply(folder.listing, files[i].secondaryFiles);
        delete files[i].secondaryFiles;
      }
      folder.listing.push (files[i])
    }
    return { "folder": folder };
  }

doc: |
  Returns Directory object with a listing set from File or File[]. Only one nested level of secondaryFiles
  is supported