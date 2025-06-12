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
    default: "#!/bin/bash\nORIGINAL_FASTQ_1=$0\nTRIMMED_FASTQ_1=$1\nTRIMMING_REPORT_1=$2\n\nORIGINAL_FASTQ_2=$3\nTRIMMED_FASTQ_2=$4\nTRIMMING_REPORT_2=$5\n\nMIN_COUNTS=$6\nTRIMMED_COUNTS=$(( `cat $TRIMMED_FASTQ_1 | wc -l` / 4 )) \n\necho \"ORIGINAL_FASTQ_1: ${ORIGINAL_FASTQ_1}\"\necho \"TRIMMED_FASTQ_1: ${TRIMMED_FASTQ_1}\"\necho \"TRIMMING_REPORT_1: ${TRIMMING_REPORT_1}\"\n\necho \"ORIGINAL_FASTQ_2: ${ORIGINAL_FASTQ_2}\"\necho \"TRIMMED_FASTQ_2: ${TRIMMED_FASTQ_2}\"\necho \"TRIMMING_REPORT_2: ${TRIMMING_REPORT_2}\"\n\necho \"TRIMMED_COUNTS: ${TRIMMED_COUNTS}\"\n\nif (( $TRIMMED_COUNTS < $MIN_COUNTS )); then\n  echo \"Bypassing adapter trimming\"\n  cp $ORIGINAL_FASTQ_1 `basename $TRIMMED_FASTQ_1`\n  cp $ORIGINAL_FASTQ_2 `basename $TRIMMED_FASTQ_2`\nelse\n  echo \"Using adapter trimming results\"\n  cp $TRIMMED_FASTQ_1 .\n  cp $TRIMMED_FASTQ_2 .\n  cp $TRIMMING_REPORT_1 .\n  cp $TRIMMING_REPORT_2 .\nfi\n"
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
baseCommand:
- bash
- -c
doc: |
  If the number of reads in the trimmed_fastq_file_1 is less then min_reads_count, tool
  will return original_fastq_file_1/2 and nulls as selected_report_file_1/2. Otherwise,
  the trimmed_fastq_file_1/2 and trimming_report_file_1/2 will be returned. Might be
  usefull in case of trimgalore removed all reads from the original_fastq_file_1/2.
label: bypass-trimgalore-pe
