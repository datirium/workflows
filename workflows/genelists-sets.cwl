cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  genelists:
    - "filter-deseq-for-heatmap.cwl"
    - "filter-diffbind-for-heatmap.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  genelist_A_name:
    type:
    - "null"
    - string
    default: null
    label: "Genelist A"
    sd:preview:
      position: 2
    doc: "genelist A alias/sample name"
    'sd:upstreamSource': "genelists/alias"
    'sd:localLabel': true

  genelist_B_name:
    type:
    - "null"
    - string
    default: null
    label: "Genelist B:"
    sd:preview:
      position: 3
    doc: "genelist B alias/sample name"
    'sd:upstreamSource': "genelists/alias"
    'sd:localLabel': true

  filtered_list_A:
    type:
    - "null"
    - File
    default: null
    format: "http://edamontology.org/format_3475"
    label: "Filtered differential genelist A"
    doc: "filtered differential genelists from DESeq or diffbind pipelines"
    'sd:upstreamSource': "genelists/filtered_file"

  filtered_list_B:
    type:
    - "null"
    - File
    default: null
    format: "http://edamontology.org/format_3475"
    label: "Filtered differential geneslist B"
    doc: "filtered differential genelists from DESeq or diffbind pipelines"
    'sd:upstreamSource': "genelists/filtered_file"

  set_operation:
    type:
    - "null"
    - type: enum
      name: "What set result are you looking for?"
      symbols:
      - Intersection
      - Union
      - Symmetric_Difference
      - Relative_Complement
    label: "Select set operation:"
    sd:preview:
      position: 4
    doc: "Set Examples where list A = {1, 2, 3, 4} and list B = {3, 4, 5, 6}:\n
       - Intersection: A ∩ B = {3, 4}; scores from list A are reported\n
       - Union: A ∪ B = {1, 2, 3, 4, 5, 6}; for overlapping elements, scores from list A are reported\n
       - Symmetric_Difference: A △ B = {1, 2, 5, 6}\n
       - Relative_Complement: A \ B = {1, 2}"
    'sd:localLabel': true


outputs:

  genelist_filtered_set:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Filtered differentially expressed genes"
    doc: "Regions of interest formatted as headerless BED file with [chrom start end name score strand]"
    outputSource: set_operator/genelist_filtered_set
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Set results'
        Title: 'Set table'

  filtering_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Filtering stdout log"
    doc: "Filtering stdout log"
    outputSource: set_operator/log_file_stdout

  filtering_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Filtering stderr log"
    doc: "Filtering stderr log"
    outputSource: set_operator/log_file_stderr


steps:

  set_operator:
    run: ../tools/genelists-sets.cwl
    in:
      filtered_list_A: filtered_list_A
      filtered_list_B: filtered_list_B
      set_operation:
        source: set_operation
        valueFrom: $(self)
    out:
      - genelist_filtered_set
      - log_file_stdout
      - log_file_stderr


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Set Operations for filtered gene lists"
label: "Set Operations for filtered gene lists"
s:alternateName: "Set Operations for filtered gene lists"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/genelists-sets.cwl
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
  # Set Operations for filtered gene lists

  This workflow takes as input 2 filtered genelists samples and performs the user-selected set operation on them.
  The output is a single filtered gene list in the same format as the input files (headerless BED file with [chrom start end name score strand]).
  The returned score value (column 5) is always derived from file A.