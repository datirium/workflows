cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: ResourceRequirement
  ramMin: 7024
  coresMin: 1
hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-kraken2:v1.0.0
inputs:
  script_command:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\n\nStdout log file for k2-download-db.cwl tool:\n\n"
      DATABASE=$0;
      printf "INPUTS:\n\n"
      printf "\$0 - $DATABASE\n\n"
      printf "EXECUTION:\n\n"
      # commands start
      if [[ "$DATABASE" == "Viral" ]]; then
        url="https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240904.tar.gz"
        db="k2_viral_20240904"
      elif [[ "$DATABASE" == "Standard" ]]; then
        url="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz"
        db="k2_standard_20240904"
      elif [[ "$DATABASE" == "Standard-16" ]]; then
        url="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20240904.tar.gz"
        db="k2_standard_16gb_20240904"
      elif [[ "$DATABASE" == "MinusB" ]]; then
        url="https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20240904.tar.gz"
        db="k2_minusb_20240904"
      elif [[ "$DATABASE" == "PlusPFP-16" ]]; then
        url="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20240904.tar.gz"
        db="k2_pluspfp_16gb_20240904"
      elif [[ "$DATABASE" == "EuPathDB46" ]]; then
        url="https://genome-idx.s3.amazonaws.com/kraken/k2_eupathdb48_20230407.tar.gz"
        db="k2_eupathdb48_20230107"
      elif [[ "$DATABASE" == "16S_Greengenes" ]]; then
        url="https://genome-idx.s3.amazonaws.com/kraken/16S_Greengenes13.5_20200326.tgz"
        db="k2_16S_Greengenes_20200326"
      elif [[ "$DATABASE" == "16S_Silva_138" ]]; then
        url="https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz"
        db="k2_16S_Silva_138_20200326"
      fi
      printf "Downloading Kraken2 $DATABASE database from: $url\n\n"
      wget $url
      mkdir ./k2db
      tar -xf *.tar.gz -C ./k2db
      printf "END OF SCRIPT\n\n"
    inputBinding:
      position: 4
  user_selection:
    type: string
    label: Name of kraken2 database to download
    inputBinding:
      position: 5
    doc: Name of kraken2 database to download and return path as output for use as an upstream input for kraken2 classify.
outputs:
  k2db:
    type: Directory
    outputBinding:
      glob: k2db
  compressed_k2db_tar:
    type: File
    outputBinding:
      glob: '*.tar.gz'
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
  Tool downloads user-specified kraken2 database from https://benlangmead.github.io/aws-indexes/k2.
  Resulting directory is used as upstream input for kraken2 classify tools.
label: k2-download-db
