cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 7024                    # equal to ~8GB
    coresMin: 1

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.3


inputs:

  script_command:
    type: string?
    default: |
      #!/bin/bash
      exec 1> error_msg.txt 2>&1
      printf "extract-fastq.cwl\n$(date)\n"
      printf "INPUTS:\n"
      printf "\$0 - $0\n"
      printf "\$1 - $1\n\n"

      shopt -s nocaseglob

      COMBINED="$0".fastq

      function extract {
        FILE=$1
        COMBINED=$2
        T=`file -b "${FILE}" | awk '{print $1}'`
        case "${T}" in
          "bzip2"|"gzip"|"Zip")
            7z e -so "${FILE}" >> "${COMBINED}"
            ;;
          "ASCII")
            cat "${FILE}" >> "${COMBINED}" || true
            ;;
          *)
            echo "Error: file type '${T}' unknown" > error_report.txt
            rm -f "${COMBINED}"
            exit 1
        esac
      } 

      for FILE in "$@"; do
          echo "Extracting:" $FILE;
          extract "${FILE}" "${COMBINED}"
      done;

    inputBinding:
      position: 5
    doc: "Bash script to extract compressed FASTQ file"

  output_prefix:
    type: string?
    inputBinding:
      position: 6
    default: "merged"
    doc: "Output prefix for extracted file"

  compressed_file:
    type:
    - File
    - type: array
      items: File
    inputBinding:
      position: 7
    doc: "Compressed or uncompressed FASTQ file(s)"


outputs:

  error_msg:
    type: File?
    outputBinding:
      glob: "error_msg.txt"

  error_report:
    type: File?
    outputBinding:
      glob: "error_report.txt"

  fastq_file:
    type: File
    outputBinding:
      glob: "*.fastq"


baseCommand: [bash, '-c']


label: "extract-fastq"
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
