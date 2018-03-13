#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_target_name = function() {
        return inputs.target_filename.split('/').slice(-1)[0];
    }

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.2

inputs:
  source_file:
    type:
      - File
    inputBinding:
      position: 5

  target_filename:
    type: string
    inputBinding:
      position: 6
      valueFrom: $(get_target_name())

outputs:
  target_file:
    type: File
    outputBinding:
      glob: $(get_target_name())


baseCommand: [cp]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "rename"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/rename.cwl
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
  Tool renames `source_file` to `target_filename`.
  Input `target_filename` shoudl be set as string. If it's a full path, only basename will be used.

s:about: |
  cp source target