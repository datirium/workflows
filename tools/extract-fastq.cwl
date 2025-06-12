cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ShellCommandRequirement
- class: ResourceRequirement
  ramMin: 7024
  coresMin: 1
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.3
inputs:
  script:
    type: string?
    default: "#!/bin/bash\nshopt -s nocaseglob\n\nCOMBINED=\"$0\".fastq\n\nfunction extract {\n  FILE=$1\n  COMBINED=$2\n  T=`file -b \"${FILE}\" | awk '{print $1}'`\n  case \"${T}\" in\n    \"bzip2\"|\"gzip\"|\"Zip\")\n      7z e -so \"${FILE}\" >> \"${COMBINED}\"\n      ;;\n    \"ASCII\")\n      cat \"${FILE}\" >> \"${COMBINED}\" || true\n      ;;\n    *)\n      echo \"Error: file type unknown\"\n      rm -f \"${COMBINED}\"\n      exit 1\n  esac\n} \n\nfor FILE in \"$@\"; do\n    echo \"Extracting:\" $FILE;\n    extract \"${FILE}\" \"${COMBINED}\"\ndone;\n"
    inputBinding:
      position: 5
    doc: |
      Bash script to extract compressed FASTQ file
  output_prefix:
    type: string?
    inputBinding:
      position: 6
    default: merged
    doc: |
      Output prefix for extracted file
  compressed_file:
    type:
    - File
    - type: array
      items: File
    inputBinding:
      position: 7
    doc: |
      Compressed or uncompressed FASTQ file(s)
outputs:
  fastq_file:
    type: File
    outputBinding:
      glob: '*'
baseCommand:
- bash
- -c
doc: |
  Tool to decompress input FASTQ file(s).
  If several FASTQ files are provided, they will be concatenated in the order that corresponds to files in input.
  Bash script's logic:
  - disable case sensitive glob check
  - check if root name of input file already include '.fastq' or '.fq' extension. If yes, set DEFAULT_EXT to "",
    otherwise use '.fastq'
  - check file type, decompress if needed
  - return 1, if file type is not recognized
  This script also works of input file doesn't have any extension at all
label: extract-fastq
