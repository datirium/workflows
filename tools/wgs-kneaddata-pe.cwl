cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 30510                                     # equal to ~32GB
    coresMin: 8


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-metaphlan:v1.0.0


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\n\nStdout log file for wgs-kneaddata-pe.cwl tool:\n\n"
      R1=$0
      R2=$1
      printf "EXECUTION:\n\n"
      # commands start
      #   rename input files for kneaddata
      cp $R1 input_for_kd.1.fastq
      cp $R2 input_for_kd.2.fastq
      mkdir kd_output
      kneaddata -t 8 --input1 input_for_kd.1.fastq --input2 input_for_kd.2.fastq -db /dockerdata/databases/human_genome_bowtie2/ --bypass-trim -o kd_output
      echo "cleaning up tmp files"
      rm input_for_kd.1.fastq
      rm input_for_kd.2.fastq
      printf "END OF SCRIPT\n\n"
    inputBinding:
        position: 4

  read1file:
    type: File
    label: "R1 fastq"
    inputBinding:
      position: 6
    doc: "FASTQ file 1 of paired end read data."

  read2file:
    type: File
    label: "R2 fastq"
    inputBinding:
      position: 7
    doc: "FASTQ file 2 of paired end read data."


outputs:

  kneaddata_cleaned_R1:
    type: File
    outputBinding:
      glob: "kd_output/*_paired_1.fastq"

  kneaddata_cleaned_R2:
    type: File
    outputBinding:
      glob: "kd_output/*_paired_2.fastq"

  kneaddata_contaminated_R1:
    type: File
    outputBinding:
      glob: "kd_output/*_paired_contam_1.fastq"

  kneaddata_contaminated_R2:
    type: File
    outputBinding:
      glob: "kd_output/*_paired_contam_2.fastq"

  kneaddata_log:
    type: File
    outputBinding:
      glob: "kd_output/*_kneaddata.log"

  stdout_log:
    type: File
    outputBinding:
      glob: "log.stdout"
    doc: |
      log for stdout

  stderr_log:
    type: File
    outputBinding:
      glob: "log.stderr"
    doc: |
      log for stderr


baseCommand: ["bash", "-c"]
stdout: log.stdout
stderr: log.stderr


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "wgs-kneaddata-pe"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/wgs-kneaddata-pe.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium LLC"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: ""
    s:streetAddress: ""
    s:telephone: ""
  s:logo: "https://avatars.githubusercontent.com/u/33202955?s=200&v=4"
  s:department:
  - class: s:Organization
    s:legalName: "Datirium LLC"
    s:department:
    - class: s:Organization
      s:legalName: "Bioinformatics"
      s:member:
      - class: s:Person
        s:name: Robert Player
        s:email: mailto:support@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X

doc: |
    Tool filters out human contamination from input read files.
    Outputs both contaminate labeled read pairs and cleaned read pairs.
