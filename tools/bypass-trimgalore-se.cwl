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
      printf "bypass-trimgalore-se.cwl\n$(date)\n"

      ORIGINAL_FASTQ=$0
      TRIMMED_FASTQ=$1
      TRIMMING_REPORT=$2
      MIN_COUNTS=$3
      TRIMMED_COUNTS=$(( `cat $TRIMMED_FASTQ | wc -l` / 4 )) 
      
      echo "ORIGINAL_FASTQ: ${ORIGINAL_FASTQ}"
      echo "TRIMMED_FASTQ: ${TRIMMED_FASTQ}"
      echo "TRIMMING_REPORT: ${TRIMMING_REPORT}"
      echo "TRIMMED_COUNTS: ${TRIMMED_COUNTS}"

      if (( $TRIMMED_COUNTS < $MIN_COUNTS )); then
        echo "Bypassing adapter trimming"
        cp $ORIGINAL_FASTQ `basename $TRIMMED_FASTQ`
      else
        echo "Using adapter trimming results"
        cp $TRIMMED_FASTQ .
        cp $TRIMMING_REPORT .
      fi
    inputBinding:
      position: 5

  original_fastq_file:
    type: File
    inputBinding:
      position: 6

  trimmed_fastq_file:
    type: File
    inputBinding:
      position: 7

  trimming_report_file:
    type: File
    inputBinding:
      position: 8

  min_reads_count:
    type: int?
    default: 100000
    inputBinding:
      position: 9


outputs:

  selected_fastq_file:
    type: File
    outputBinding:
      glob: $(inputs.trimmed_fastq_file.basename)
      
  selected_report_file:
    type: File?
    outputBinding:
      glob: $(inputs.trimming_report_file.basename)


baseCommand: ["bash", "-c"]


label: "bypass-trimgalore-se"
doc: |
  If the number of reads in the trimmed_fastq_file is less then min_reads_count, tool
  will return original_fastq_file and null as selected_report_file. Otherwise, the
  trimmed_fastq_file and trimming_report_file will be returned. Might be usefull in
  case of trimgalore removed all reads from the original_fastq_file