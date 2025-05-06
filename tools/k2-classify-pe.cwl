cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 244130                                     # equal to ~264GB
    coresMin: 20


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-kraken2:v1.0.0


inputs:

  k2db:
    type: Directory
    inputBinding:
      prefix: "-d"
    doc: "Name of kraken2 database to use for kraken2 classify."

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

  k2_classified_R1:
    type: File
    outputBinding:
      glob: "k2_classified_reads_1.fastq"

  k2_classified_R2:
    type: File
    outputBinding:
      glob: "k2_classified_reads_2.fastq"

  k2_unclassified_R1:
    type: File
    outputBinding:
      glob: "k2_unclassified_reads_1.fastq"

  k2_unclassified_R2:
    type: File
    outputBinding:
      glob: "k2_unclassified_reads_2.fastq"

  k2_output:
    type: File
    outputBinding:
      glob: "k2.output"

  k2_report:
    type: File
    outputBinding:
      glob: "k2.report"

  k2_report_tsv:
    type: File
    outputBinding:
      glob: "k2_report.tsv"

  k2_stderr:
    type: File
    outputBinding:
      glob: "parsed.stderr"

  krona_html:
    type: File?
    outputBinding:
      glob: "krona.html"


baseCommand: ["/usr/local/bin/run_kraken2pe.sh"]
stdout: error_report.txt
stderr: error_msg.txt


label: "k2-classify-pe"
doc: |
  Taxonomically classifies paired-end sequencing reads in FASTQ format, that have been
  adapter trimmed with trimgalore, using Kraken2 with a user-selected pre-built database
