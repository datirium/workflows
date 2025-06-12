cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4
inputs:
  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\nLog file for samtools-clean-headers.cwl tool:\n\n"
      bam=$0;
      samtools view -H $bam > clean.sam
      samtools view -S $bam | awk -F'\t' 'BEGIN{OFS="\t"};{sub(/ .*/,"",$1); print($0)}' >> clean.sam
      samtools view -hb clean.sam > clean.bam
      # delete any intermediate files from system
      rm clean.sam
    inputBinding:
      position: 4
  bam_file:
    type: File
    label: Input BAM file
    inputBinding:
      position: 5
    doc: Input BAM file, does not have to be coordinates sorted
outputs:
  preseq_bam:
    type: File?
    outputBinding:
      glob: clean.bam
  log_file_stdout:
    type: File
    outputBinding:
      glob: samtools-clean-headers.log.stdout
    doc: |
      log for stdout
  log_file_stderr:
    type: File
    outputBinding:
      glob: samtools-clean-headers.log.stderr
    doc: |
      log for stderr
baseCommand:
- bash
- -c
stdout: samtools-clean-headers.log.stdout
stderr: samtools-clean-headers.log.stderr
doc: |
  Tool processes BAM file, to remove everything in the header (col1) after a single
  space is found (including the space). This cleaned bam file is then used as input
  for the preseq step (tools/preseq-lc-extrap.cwl).
label: samtools-clean-headers
