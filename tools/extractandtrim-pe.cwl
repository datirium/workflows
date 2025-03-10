cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 15250                                     # equal to ~16GB
    coresMin: 10


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-extractandtrim:dev


inputs:

  read1file:
    type: File
    inputBinding:
      prefix: "-a"
    doc: "FASTQ file 1 of paired end read data."

  read2file:
    type: File
    inputBinding:
      prefix: "-b"
    doc: "FASTQ file 2 of paired end read data."


outputs:

  error_msg:
    type: File?
    outputBinding:
      glob: "error_msg.txt"

  error_report:
    type: File?
    outputBinding:
      glob: "error_report.txt"

  trimmed_fastq_r1:
    type: File
    outputBinding:
      glob: "extracted_combined_R1_val_1.fq"

  trimmed_fastq_r2:
    type: File
    outputBinding:
      glob: "extracted_combined_R2_val_2.fq"

  trimmed_fastq_r1_report:
    type: File
    outputBinding:
      glob: "extracted_combined_R1.fastq_trimming_report.txt"

  trimmed_fastq_r2_report:
    type: File
    outputBinding:
      glob: "extracted_combined_R2.fastq_trimming_report.txt"

  statistics_file_r1:
    type: File
    outputBinding:
      glob: "fastx_quality_stats_r1.tsv"

  statistics_file_r2:
    type: File
    outputBinding:
      glob: "fastx_quality_stats_r2.tsv"


baseCommand: ["/usr/local/bin/run_extractandtrimpe.sh"]
stdout: error_report.txt
stderr: error_msg.txt


label: "extractandtrim-pe"
doc: |
  Extracts compressed fastq paired-end read files and concatentates each pair if multiple files.
  Runs each concatenated read pair file through trim galore, then generates summarized read statistics for the trimmed reads.
