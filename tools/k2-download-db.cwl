cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 7024                                     # equal to ~8GB
    coresMin: 1


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-kraken2:v1.0.0


inputs:

  user_selection:
    type: string
    inputBinding:
      prefix: "-d"
    doc: "kraken2 database name to download (acceptable strings are: Viral, Standard, Standard-16, MinusB, PlusPFP-16, EuPathDB46, 16S_Greengenes, 16S_Silva_138)"


outputs:

  error_msg:
    type: File?
    outputBinding:
      glob: "error_msg.txt"

  error_report:
    type: File?
    outputBinding:
      glob: "error_report.txt"

  k2db:
    type: Directory
    outputBinding:
      glob: "k2db"

  compressed_k2db_tarball:
    type: File
    outputBinding:
      glob: "*.tar.gz"


baseCommand: ["/usr/local/bin/run_kraken2download.sh"]
stdout: error_report.txt
stderr: error_msg.txt


label: "k2-download-db"
doc: |
    Tool downloads user-specified kraken2 database from https://benlangmead.github.io/aws-indexes/k2.
    Resulting directory is used as upstream input for kraken2 classify tools.
