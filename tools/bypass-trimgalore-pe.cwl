cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: ResourceRequirement
    ramMin: 7024
    coresMin: 1


hints:
- class: DockerRequirement
  dockerPull: scidap/scidap:v0.0.4


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      exec 1> error_msg.txt 2>&1
      printf "bypass-trimgalore-pe.cwl\n$(date)\n"

      ORIGINAL_FASTQ_1=$0
      TRIMMED_FASTQ_1=$1
      TRIMMING_REPORT_1=$2

      ORIGINAL_FASTQ_2=$3
      TRIMMED_FASTQ_2=$4
      TRIMMING_REPORT_2=$5

      MIN_COUNTS=$6
      TRIMMED_COUNTS=$(( `cat $TRIMMED_FASTQ_1 | wc -l` / 4 )) 
      
      echo "ORIGINAL_FASTQ_1: ${ORIGINAL_FASTQ_1}"
      echo "TRIMMED_FASTQ_1: ${TRIMMED_FASTQ_1}"
      echo "TRIMMING_REPORT_1: ${TRIMMING_REPORT_1}"

      echo "ORIGINAL_FASTQ_2: ${ORIGINAL_FASTQ_2}"
      echo "TRIMMED_FASTQ_2: ${TRIMMED_FASTQ_2}"
      echo "TRIMMING_REPORT_2: ${TRIMMING_REPORT_2}"

      echo "TRIMMED_COUNTS: ${TRIMMED_COUNTS}"

      if (( $TRIMMED_COUNTS < $MIN_COUNTS )); then
        echo "Bypassing adapter trimming"
        cp $ORIGINAL_FASTQ_1 `basename $TRIMMED_FASTQ_1`
        cp $ORIGINAL_FASTQ_2 `basename $TRIMMED_FASTQ_2`
      else
        echo "Using adapter trimming results"
        cp $TRIMMED_FASTQ_1 .
        cp $TRIMMED_FASTQ_2 .
        cp $TRIMMING_REPORT_1 .
        cp $TRIMMING_REPORT_2 .
      fi
    inputBinding:
      position: 5

  original_fastq_file_1:
    type: File
    inputBinding:
      position: 6

  trimmed_fastq_file_1:
    type: File
    inputBinding:
      position: 7

  trimming_report_file_1:
    type: File
    inputBinding:
      position: 8

  original_fastq_file_2:
    type: File
    inputBinding:
      position: 9

  trimmed_fastq_file_2:
    type: File
    inputBinding:
      position: 10

  trimming_report_file_2:
    type: File
    inputBinding:
      position: 11

  min_reads_count:
    type: int?
    default: 100000
    inputBinding:
      position: 12


outputs:

  error_msg:
    type: File?
    outputBinding:
      glob: "error_msg.txt"

  error_report:
    type: File?
    outputBinding:
      glob: "error_report.txt"

  selected_fastq_file_1:
    type: File
    outputBinding:
      glob: $(inputs.trimmed_fastq_file_1.basename)
      
  selected_report_file_1:
    type: File?
    outputBinding:
      glob: $(inputs.trimming_report_file_1.basename)

  selected_fastq_file_2:
    type: File
    outputBinding:
      glob: $(inputs.trimmed_fastq_file_2.basename)
      
  selected_report_file_2:
    type: File?
    outputBinding:
      glob: $(inputs.trimming_report_file_2.basename)


baseCommand: ["bash", "-c"]


label: "bypass-trimgalore-pe"
doc: |
  If the number of reads in the trimmed_fastq_file_1 is less then min_reads_count, tool
  will return original_fastq_file_1/2 and nulls as selected_report_file_1/2. Otherwise,
  the trimmed_fastq_file_1/2 and trimming_report_file_1/2 will be returned. Might be
  usefull in case of trimgalore removed all reads from the original_fastq_file_1/2.