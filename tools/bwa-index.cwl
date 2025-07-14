cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ShellCommandRequirement
hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-gatk4:v1.0.0
inputs:
  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\n\nStdout log file for bwa-index.cwl tool:\n\n"
      FASTA=$0;
      printf "INPUTS:\n\n"
      printf "\$0 - $FASTA\n\n"
      printf "EXECUTION:\n\n"
      # commands start
      mkdir ./bwa-index
      cd bwa-index
      cp $FASTA genome.fasta
      bwa index genome.fasta
      # ref file need fasta `.fai` index for HaplotypeCaller
      samtools faidx genome.fasta
      samtools dict genome.fasta > genome.dict
    inputBinding:
      position: 1
  ref_genome_fasta:
    type: File
    inputBinding:
      position: 2
    doc: |
      Reference genome FASTA to be indexed with 'bwa index'
outputs:
  bwa_index:
    type: Directory
    outputBinding:
      glob: bwa-index
  log_file_stdout:
    type: File
    outputBinding:
      glob: log.stdout
    doc: |
      log for stdout
  log_file_stderr:
    type: File
    outputBinding:
      glob: log.stderr
    doc: |
      log for stderr
baseCommand:
- bash
- -c
stdout: log.stdout
stderr: log.stderr
doc: |
  Tool indexes a reference genome fasta using bwa index. Currently these indexes is used in the
  germline variant calling workflow.
label: bwa-index
