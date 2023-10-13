cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-kallisto:stable


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
    label: "threads"
    inputBinding:
      position: 3
    doc: "Number of threads for steps that support multithreading."


outputs:

  kallisto_index:
    type: File
    outputBinding:
      glob: kallisto-index.kdx

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

s:name: "kallisto-index"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/kallisto-index.cwl
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
    Tool indexes a reference genome fasta using `kallisto index`.
