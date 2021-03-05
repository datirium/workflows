cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  sample_to_filter:
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

  feature_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "ChIP/ATAC experiment"
    doc: "Called peaks file in TSV format with the nearest genes assigned"
    'sd:upstreamSource': "sample_to_filter/iaintersect_result"
    'sd:localLabel': true

  alignment_file:
    type: File
    secondaryFiles:
    - .bai
    format: "http://edamontology.org/format_2572"
    label: "ChIP/ATAC experiment"
    doc: "Coordinate sorted BAM and BAI index files"
    'sd:upstreamSource': "sample_to_filter/bambai_pair"

  sql_query:
    type: string
    label: "Filtering parameters"
    doc: "Filtering parameters (WHERE parameters for SQL query)"
    'sd:filtering':
      params:
        columns: ["refseq_id", "gene_id", "txStart", "txEnd", "strand", "chrom", "start", "end", "length", "abssummit", "pileup", "log10p", "foldenrich", "log10q", "region"]
        types:   ["string", "string", "number", "number","string", "string","number", "number", "number", "number", "number", "number", "number", "number", "string"]

  header:
    type: boolean?
    default: false
    label: "Include header line"
    doc: "Print header line in the output file"
    'sd:layout':
      advanced: true

  columns:
    type:
    - "null"
    - string[]
    default: ["chrom", "start", "end", "chrom || '-' || start || '-' || end AS name"]
    label: "Columns to print"
    doc: |
      List of columns to print (SELECT parameters for SQL query).
      Need to have format [chrom start end name]. No header.
      4th columns should be unique, so we use a combination of chrom-start-end.
    'sd:layout':
      advanced: true


outputs:

  filtered_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Filtered called peaks with the nearest genes assigned"
    doc: "Regions of interest formatted as headerless BED file with [chrom start end name]"
    outputSource: feature_select/filtered_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Filtering results'
        Title: 'Filtered table'

  bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Not changed coordinate sorted BAM and BAI index files"
    doc: "Not changed coordinate sorted BAM and BAI index files"
    outputSource: alignment_file

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
      feature_file: feature_file
      sql_query: sql_query
      columns:
        source: columns
        valueFrom: $("DISTINCT " + self.join(", "))   # multiple peaks can have the same coordinates but different abssummit, so we need to use DISTINCT
      header: header
    out:
    - filtered_file
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Filter ChIP/ATAC peaks for Tag Density Profile or Motif Enrichment analyses"
label: "Filter ChIP/ATAC peaks for Tag Density Profile or Motif Enrichment analyses"
s:alternateName: "Filter ChIP/ATAC peaks for Tag Density Profile or Motif Enrichment analyses"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/filter-peaks-for-heatmap.cwl
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
  Filters ChIP/ATAC peaks with the neatest genes assigned for Tag Density Profile or Motif Enrichment analyses
  ============================================================================================================

  Tool filters output from any ChIP/ATAC pipeline to create a file with regions of interest for Tag Density
  Profile or Motif Enrichment analyses. Peaks with duplicated coordinates are discarded.
