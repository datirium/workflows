cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  sample_to_filter:
    - "deseq.cwl"
    - "deseq-lrt.cwl"
    - "chipseq-se.cwl"
    - "chipseq-pe.cwl"
    - "trim-chipseq-se.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-se.cwl"
    - "trim-atacseq-pe.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  diff_expr_file:
    type: File?
    default: null
    format: "http://edamontology.org/format_3475"
    label: "Input file to filter"
    doc: "Output from any DESeq like pipeline"
    'sd:upstreamSource': "sample_to_filter/diff_expr_file"
    'sd:localLabel': true

  iaintersect_result:
    type: File?
    default: null
    format: "http://edamontology.org/format_3475"
    label: "Input file to filter"
    doc: "Output from any ChiP-Seq like pipeline"
    'sd:upstreamSource': "sample_to_filter/iaintersect_result"

  sql_query:
    type: string
    label: "Filtering parameters"
    doc: "Filtering parameters to be appended after 'SELECT * FROM _ WHERE' statement"
    'sd:filtering':
    - diff_expr_file
    - iaintersect_result

  header:
    type: boolean?
    default: false
    label: "Print header line in the output file"
    doc: "Print header line in the output file"
    'sd:layout':
      advanced: true

  columns:
    type: string?
    default: "*"
    label: "Comma-separated list of column names to print (or WHERE parameters for SQL query)"
    doc: "Comma-separated list of column names to print. Default: all"
    'sd:layout':
      advanced: true


outputs:

  filtered_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Filtered TSV file"
    doc: "Filtered by provided SQL query TSV file"
    outputSource: feature_select/filtered_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Filtering results'
        Title: 'Filtered table'

  filtering_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Filtering stdout log"
    doc: "Filtering stdout log"
    outputSource: feature_select/stdout_log

  filtering_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Filtering stderr log"
    doc: "Filtering stderr log"
    outputSource: feature_select/stderr_log


steps:

  feature_select:
    run: ../tools/feature-select-sql.cwl
    in:
      feature_file:
        source: [diff_expr_file, iaintersect_result]
        valueFrom: $(self[0] || self[1])
      sql_query: sql_query
      columns: columns
      header: header
    out:
    - filtered_file
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Feature select - filters TSV/CSV files based on the provided SQL query parameters"
label: "Feature select - filters TSV/CSV files based on the provided SQL query parameters"
s:alternateName: "Feature select - filters TSV/CSV files based on the provided SQL query parameters"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/feature-select-sql.cwl
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
  Feature select - filters TSV/CSV files based on the provided SQL query parameters
  =================================================================================

  Tool filters input TSV/CSV file based on the provided SQL query.
  Value set in sql_query parameter will be appended to the
  "SELECT * FROM raw_data WHERE". Column names in sql_query are case
  insensitive. Format of the input files is identified based on file
  extension
  *.csv - CSV
  *.tsv - TSV
  Otherwise used CSV by default. Output is always saved in TSV format.
