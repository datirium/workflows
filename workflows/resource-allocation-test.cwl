cwlVersion: v1.0
class: Workflow


inputs:

  alias:
    type: string
    label: "Analysis name"
    sd:preview:
      position: 1


  delay_1:
    type: int?
    default: 60
    label: "Delay 1, sec"
    doc: "Delay 1, sec"

  threads_1:
    type: int?
    default: 1
    label: "CPUs 1"
    doc: "CPUs 1"

  memory_1:
    type: int?
    default: 15259
    doc: "Memory 1"
    label: "Memory 1"


  delay_2:
    type: int?
    default: 120
    label: "Delay 2, sec"
    doc: "Delay 2, sec"

  threads_2:
    type: int?
    default: 2
    label: "CPUs 2"
    doc: "CPUs 2"

  memory_2:
    type: int?
    default: 30518
    doc: "Memory 2"
    label: "Memory 2"


  delay_3:
    type: int?
    default: 180
    label: "Delay 3, sec"
    doc: "Delay 3, sec"

  threads_3:
    type: int?
    default: 3
    label: "CPUs 3"
    doc: "CPUs 3"

  memory_3:
    type: int?
    default: 61036
    doc: "Memory 3"
    label: "Memory 3"


outputs:

  dummy_1:
    type: File?
    outputSource: sleep_1/dummy

  dummy_2:
    type: File?
    outputSource: sleep_2/dummy

  dummy_3:
    type: File?
    outputSource: sleep_3/dummy


steps:

  sleep_1:
    run: ../tools/sleep.cwl
    in:
      delay: delay_1
      threads: threads_1
      memory: memory_1
    out: [dummy]

  sleep_2:
    run: ../tools/sleep.cwl
    in:
      delay: delay_2
      threads: threads_2
      memory: memory_2
    out: [dummy]

  sleep_3:
    run: ../tools/sleep.cwl
    in:
      delay: delay_3
      threads: threads_3
      memory: memory_3
    out: [dummy]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Resource Allocation Test"
label: "Resource Allocation Test"
s:alternateName: "Resource Allocation Test"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/resource-allocation-test.cwl
s:codeRepository: https://github.com/datirium/workflows
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
  Resource Allocation Test