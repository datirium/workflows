cwlVersion: v1.0
class: Workflow


requirements:
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement


inputs:

  isoforms_file:
    type: File[]
    doc: "Isoforms CSV file(s)"


outputs:

  genes_file:
    type: File[]
    outputSource: group_isoforms/genes_file
    doc: "Gene expression file(s) in TSV format"

  common_tss_file:
    type: File[]
    outputSource: group_isoforms/common_tss_file
    doc: "Common tss expression file(s) in TSV format"


steps:

  group_isoforms:
    run: ../tools/group-isoforms.cwl
    in:
      isoforms_file: isoforms_file
    scatter:
      - isoforms_file
    out:
      - genes_file
      - common_tss_file


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "group-isoforms-batch"
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
  Workflow runs group-isoforms.cwl tool using scatter for isoforms_file input. genes_filename and common_tss_filename inputs are ignored.

s:about: |
  Workflow runs group-isoforms.cwl tool using scatter for isoforms_file input. genes_filename and common_tss_filename inputs are ignored.