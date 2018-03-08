#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool


requirements:
- $import: ./metadata/envvar-global.yml
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.2


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      FILE=$0
      T=`file -b "${FILE}" | awk '{print $1}'`
      case "${T}" in
        "bzip2"|"gzip"|"Zip")
          7z e "${FILE}"
          ;;
        "ASCII")
          cp "${FILE}" .
          ;;
        *)
          echo "Error: file type unknown"
          exit 1
      esac
    inputBinding:
      position: 5
    doc: |
      Bash script to extract compressed FASTQ file

  compressed_file:
    type: File
    inputBinding:
      position: 6
    doc: |
      Compressed or uncompressed FASTQ file


outputs:

  fastq_file:
    type: File
    outputBinding:
      glob: "*"

baseCommand: [bash, '-c']


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "extract-fastq"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/extract-fastq.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Tool to decompress imput FASTQ file. If file is not compressed, return original FASTQ file

s:about: |
  Tool to decompress imput FASTQ file. If file is not compressed, return original FASTQ file