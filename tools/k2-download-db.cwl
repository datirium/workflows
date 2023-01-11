cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-kraken2:dev


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\nLog file for k2-download-db.cwl tool:\n\nCommands:\n"
      DATABASE=$0;
      if [[ "$DATABASE" == "Viral" ]]; then
        url="https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz"
        db="k2_viral_20221209"
      elif [[ "$DATABASE" == "MinusB" ]]; then
        url="https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20221209.tar.gz"
        db="k2_minusb_20221209"
      elif [[ "$DATABASE" == "PlusPFP-16" ]]; then
        url="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20221209.tar.gz"
        db="k2_pluspfp_16gb_20221209"
      elif [[ "$DATABASE" == "EuPathDB46" ]]; then
        url="https://genome-idx.s3.amazonaws.com/kraken/k2_eupathdb48_20201113.tar.gz"
        db="k2_eupathdb48_20201113"
      fi
      cat << EOF
        wget $url
        mkdir ./k2db
        tar -xf *.tar.gz -C ./k2db
        rm *.tar.gz
      EOF
      printf "Execution:\n"
      wget $url
      mkdir ./k2db
      tar -xf *.tar.gz -C ./k2db
      rm *.tar.gz
    inputBinding:
        position: 4

  user_selection:
    type: string
    label: "Name of kraken2 database to download"
    inputBinding:
      position: 5
    doc: "Name of kraken2 database to download and return path as output for use as an upstream input for kraken2 classify."


outputs:

  k2db:
    type: Directory
    outputBinding:
      glob: "k2db"

  log_file_stdout:
    type: File
    outputBinding:
      glob: "log.stdout"
    doc: |
      log for stdout

  log_file_stderr:
    type: File
    outputBinding:
      glob: "log.stderr"
    doc: |
      log for stderr


baseCommand: ["bash", "-c"]
stdout: 'log.stdout'
stderr: 'log.stderr'


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "k2-download-db"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/k2-download-db.cwl
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
    Tool downloads specified kraken2 database from https://benlangmead.github.io/aws-indexes/k2.
    Resulting directory is used as upstream input to kraken2 classify tools.
