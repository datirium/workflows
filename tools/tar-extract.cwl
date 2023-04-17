cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: ubuntu:20.04


inputs:

  file_to_extract:
    type: File
    inputBinding:
      position: 1
    doc: "File to extract"


outputs:

  extracted_folder:
    type: Directory
    outputBinding:
      glob: "*"
    doc: "Extracted folder"


baseCommand: ["tar", "xzf"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "TAR extract"
s:name: "TAR extract"
s:alternateName: "Extracts the content of TAR file into a folder"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/tar-extract.cwl
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
  TAR extract
  ===============================================
  
  Extracts the content of TAR file into a folder.