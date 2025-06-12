cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.3
inputs:
  script:
    type: string?
    default: "#!/bin/bash\nCOMBINED=\"$0\"\nfunction extract {\n  FILE=$1\n  COMBINED=$2\n  T=`file -b \"${FILE}\" | awk '{print $1}'`\n  case \"${T}\" in\n    \"bzip2\"|\"gzip\"|\"Zip\")\n      7z e -so \"${FILE}\" >> \"${COMBINED}\"\n      ;;\n    \"ASCII\")\n      cat \"${FILE}\" >> \"${COMBINED}\" || true\n      ;;\n    *)\n      echo \"Error: file type unknown\"\n      rm -f \"${COMBINED}\"\n      exit 1\n  esac\n} \nfor FILE in \"$@\"; do\n    echo \"Extracting:\" $FILE;\n    extract \"${FILE}\" \"${COMBINED}\"\ndone;\n"
    inputBinding:
      position: 5
    doc: |
      Bash script to extract compressed file(s)
  output_filename:
    type: string
    inputBinding:
      position: 6
    doc: |
      Output filename for extracted and optionally
      merged file(s)
  file_to_extract:
    type:
    - File
    - type: array
      items: File
    inputBinding:
      position: 7
    doc: |
      Compressed file(s) to extract
outputs:
  extracted_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
baseCommand:
- bash
- -c
doc: |
  Tool to decompress input file(s).
  If several files are provided, they will be concatenated in
  the order that corresponds to files in input.
  Bash script's logic:
  - check file type, decompress if needed, otherwise just cat
    the content of the file
  - return 1, if file type is not recognized
label: extract-fastq
