cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ShellCommandRequirement
hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-kallisto:v1.0.0
inputs:
  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\n\nStdout log file for kallisto-index.cwl tool:\n\n"
      FASTA=$0; THREADS=$1;
      printf "INPUTS:\n\n"
      printf "\$0 - $FASTA\n\n"
      printf "\$1 - $THREADS\n\n"
      # commands start
      #   detect and decompress
      if [[ $(basename $FASTA | sed 's/.*\.//') == "gz" || $(basename $FASTA | sed 's/.*\.//') == "fasta" || $(basename $FASTA | sed 's/.*\.//') == "fa" || $(basename $FASTA | sed 's/.*\.//') == "fna" ]]; then
        kallisto index -t $THREADS -i kallisto-index.kdx $FASTA   # ref fasta can be compressed with gzip or not
      else
        echo "reference file must end in '.gz', '.fasta', '.fa', or '.fna'"; exit
      fi
    inputBinding:
      position: 1
  ref_genome_fasta:
    type: File
    inputBinding:
      position: 2
    doc: |
      Reference genome FASTA to be indexed with 'kallisto index'
  threads:
    type: int
    label: threads
    inputBinding:
      position: 3
    doc: Number of threads for steps that support multithreading.
outputs:
  kallisto_index:
    type: File
    outputBinding:
      glob: kallisto-index.kdx
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
  Tool indexes a reference genome fasta using `kallisto index`.
label: kallisto-index
