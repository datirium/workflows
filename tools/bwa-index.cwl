cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-gatk4:dev


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

s:name: "bwa-index"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/bwa-index.cwl
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
    Tool indexes a reference genome fasta using bwa index. Currently these indexes is used in the
    germline variant calling workflow.
