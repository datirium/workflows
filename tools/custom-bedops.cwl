cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedops:v2.4.34


inputs:

  script:
    type: string?
    default: |
      cat "$0" > `basename $0`
    inputBinding:
      position: 1

  input_file:
    type:
      - File
      - File[]
    inputBinding:
      position: 2

  param:
    type:
    - string?
    - string[]
    inputBinding:
      position: 3


outputs:

  output_file:
    type: File
    outputBinding:
      glob: "*"


baseCommand: [bash, '-c']


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "custom-bedops"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/custom-bedops.cwl
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
        s:email: mailto:michael.kotliar@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Tool to run custom script set as `script`
  input with arguments from `param`. Based
  on bedops Dockerfile.

s:about: |
  Custom bash script runner