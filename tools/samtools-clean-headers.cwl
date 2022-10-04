cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\nLog file for samtools-clean-headers.cwl tool:\n\n"
      bam=$0;
      samtools view -H $bam > clean.sam
      samtools view -S $bam | awk -F'\t' 'BEGIN{OFS="\t"};{sub(/ .*/,"",$1); print($0)}' >> clean.sam
      samtools view -hb clean.sam > clean.bam
      # delete any intermediate files from system
      rm clean.sam
    inputBinding:
        position: 4

  bam_file:
    type: File
    label: "Input BAM file"
    inputBinding:
      position: 5
    doc: "Input BAM file, does not have to be coordinates sorted"


outputs:

  preseq_bam:
    type: File
    outputBinding:
      glob: "clean.bam"

  log_file_stdout:
    type: File
    outputBinding:
      glob: "samtools-clean-headers.log.stdout"
    doc: |
      log for stdout

  log_file_stderr:
    type: File
    outputBinding:
      glob: "samtools-clean-headers.log.stderr"
    doc: |
      log for stderr


baseCommand: ["bash", "-c"]
stdout: 'samtools-clean-headers.log.stdout'
stderr: 'samtools-clean-headers.log.stderr'


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "samtools-clean-headers"
s:downloadUrl: https://github.com/datirium/workflows/tree/master/workflows/tools/samtools-clean-headers.cwl
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
        s:email: mailto:robert.player@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X

doc: |
    Tool processes BAM file, to remove everything in the header (col1) after a single
    space is found (including the space). This cleaned bam file is then used as input
    for the preseq step (tools/preseq-lc-extrap.cwl).
