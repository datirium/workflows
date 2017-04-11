#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  - $import: ./metadata/envvar-global.yml
  - class: InlineJavascriptRequirement
    expressionLib:
    - var default_output_filename = function() {
          return inputs.input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+".sorted.bed";
      };

hints:
- class: DockerRequirement
  dockerPull: ubuntu:15.10

inputs:
  input:
    type: File
    inputBinding:
      position: 4

  key:
    type:
      type: array
      items: string
      inputBinding:
        prefix: "-k"
    inputBinding:
      position: 1
    doc: |
      -k, --key=POS1[,POS2]
      start a key at POS1, end it at POS2 (origin 1)

outputs:
  sorted:
    type: File
    doc: "The sorted file"
    outputBinding:
      glob: $(default_output_filename())

stdout: $(default_output_filename())

baseCommand: ["sort"]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "linux-sort"
s:downloadUrl: https://raw.githubusercontent.com/SciDAP/workflows/master/tools/linux-sort.cwl
s:codeRepository: https://github.com/SciDAP/workflows
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
        s:name: Andrey Kartashov
        s:email: mailto:Andrey.Kartashov@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0001-9102-5681

doc: |
  Tool is used to run sort command with input data
