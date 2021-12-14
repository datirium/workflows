class: Workflow
cwlVersion: v1.1


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  mixcr_sample:
  - "mixcr-pe.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  vdj_file:
    type: File[]
    label: "MiXCR experiments"
    doc: "Array of VDJ files from MiXCR pipeline"
    'sd:upstreamSource': "mixcr_sample/vdj_file"
    'sd:localLabel': true

  vdj_name:
    type: string[]
    label: "MiXCR experiments"
    doc: "Names for input VDJ files. Order corresponds to the vdj_file"
    'sd:upstreamSource': "mixcr_sample/alias"

  intersect_type:
    type:
    - "null"
    - type: enum
      symbols:
      - "strict"
      - "nt"
      - "ntV"
      - "ntVJ"
      - "aa"
      - "aaV"
      - "aaVJ"
      - "aa!nt"
    default: "aa"
    label: "Clonotype features to compare when checking if two clonotypes match"
    doc: |
      Specifies which clonotype features (CDR3 sequence, V/J segments, hypermutations)
      will be compared when checking if two clonotypes match.

  min_times_detected:
    type: int?
    default: 2
    label: "Minimal number of samples where clonotype should be detected"
    doc: |
      Minimal number of samples where clonotype should be detected.


outputs:

  combined_vdj_file:
    type: File
    label: "Combined clonotypes file in VDJTools format"
    doc: "Combined clonotypes file in VDJTools format"
    outputSource: vdjtools_join_samples/combined_vdj_file

  vdj_summary_file:
    type: File
    label: "Summary file from joning multiple VDJ files"
    doc: "Summary file from joning multiple VDJ files"
    outputSource: vdjtools_join_samples/summary_file

  vdj_metadata_file:
    type: File
    label: "Metadata file from joning multiple VDJ files"
    doc: "Metadata file from joning multiple VDJ files"
    outputSource: vdjtools_join_samples/metadata_file

  vdj_venn_diag_plot:
    type: File
    label: "Venn diagram from joning multiple VDJ files"
    doc: "Venn diagram from joning multiple VDJ files"
    outputSource: vdjtools_join_samples/venn_diag_plot


steps:

  vdjtools_join_samples:
    run: ../tools/vdjtools-join-samples.cwl
    in:
      vdj_file: vdj_file
      vdj_name: vdj_name
      intersect_type: intersect_type
      min_times_detected: min_times_detected
    out:
    - combined_vdj_file
    - summary_file
    - metadata_file
    - venn_diag_plot


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Clonotype Abundance Analysis"
label: "Clonotype Abundance Analysis"
s:alternateName: "Joins several clonotype tables together to form a joined clonotype abundance table"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/vdjtools-join.cwl
s:codeRepository: https://github.com/Barski-lab/workflows-datirium
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
  Clonotype Abundance Analysis
  ============================
  
  Joins several clonotype tables together to form a joined
  clonotype abundance table