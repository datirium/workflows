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
    default: "#!/bin/bash\nORIGINAL_FASTQ=$0\nTRIMMED_FASTQ=$1\nTRIMMING_REPORT=$2\nMIN_COUNTS=$3\nTRIMMED_COUNTS=$(( `cat $TRIMMED_FASTQ | wc -l` / 4 )) \n\necho \"ORIGINAL_FASTQ: ${ORIGINAL_FASTQ}\"\necho \"TRIMMED_FASTQ: ${TRIMMED_FASTQ}\"\necho \"TRIMMING_REPORT: ${TRIMMING_REPORT}\"\necho \"TRIMMED_COUNTS: ${TRIMMED_COUNTS}\"\n\nif (( $TRIMMED_COUNTS < $MIN_COUNTS )); then\n  echo \"Bypassing adapter trimming\"\n  cp $ORIGINAL_FASTQ `basename $TRIMMED_FASTQ`\nelse\n  echo \"Using adapter trimming results\"\n  cp $TRIMMED_FASTQ .\n  cp $TRIMMING_REPORT .\nfi\n"
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
baseCommand:
- bash
- -c
doc: |
  If the number of reads in the trimmed_fastq_file is less then min_reads_count, tool
  will return original_fastq_file and null as selected_report_file. Otherwise, the
  trimmed_fastq_file and trimming_report_file will be returned. Might be usefull in
  case of trimgalore removed all reads from the original_fastq_file
label: bypass-trimgalore-se
