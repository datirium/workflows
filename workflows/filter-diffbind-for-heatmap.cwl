cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  sample_to_filter:
  - "diffbind.cwl"
  - "https://github.com/datirium/workflows/workflows/diffbind.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  feature_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Differential Binding Analysis experiment"
    doc: "Differential binding analysis results exported as TSV"
    'sd:upstreamSource': "sample_to_filter/diffbind_report_file"
    'sd:localLabel': true

  sql_query:
    type: string
    label: "Filtering parameters"
    doc: "Filtering parameters (WHERE parameters for SQL query)"
    'sd:filtering':
      params:
        columns: ["Refseq_id", "Gene_id", "txStart", "txEnd", "Strand", "Region", "Chr", "Start", "End", "Conc", "Conc1", "Conc2", "Fold", "p-value", "FDR", "Called1", "Called2"]
        types:   ["string", "string", "number", "number", "string", "string", "string", "number", "number", "number", "number", "number", "number", "number", "number","number", "number"]
                 
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
    default: ["Chr AS chrom", "txStart AS start", "txEnd AS end", "Gene_id AS name", "Fold AS score", "Strand AS strand"]
    label: "Columns to print"
    doc: |
      List of columns to print (SELECT parameters for SQL query).
      Need to have format [chrom start end name score strand]. No header.
      4th columns should be unique, so we use GeneId for that.
      5th columns will be ignored by Tag Density pipeline, so we use log2FoldChange.
    'sd:layout':
      advanced: true


outputs:

  filtered_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Filtered regions"
    doc: "Filtered regions of interest by default formatted as headerless BED file with [Chr Start End]"
    outputSource: feature_select/filtered_file

  filtered_file_w_header:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Filtered differentially expressed genes"
    doc: "Regions of interest formatted as headered BED file with [chrom start end name score strand]"
    outputSource: add_header/filtered_file_with_header
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

  rename_header:
    run: ../tools/custom-bash.cwl
    in:
      input_file: feature_file
      script:
        default: |
          echo "Replacing header to include Conc1 and Conc2 instead of Conc_[group1] and Conc_[group2]"
          cat "$0" | grep -v "Refseq_id" | cut -f 1-17 > headerless_report.tsv
          echo -e "Refseq_id\tGene_id\ttxStart\ttxEnd\tStrand\tRegion\tChr\tStart\tEnd\tConc\tConc1\tConc2\tFold\tp-value\tFDR\tCalled1\tCalled2" > `basename $0`
          cat headerless_report.tsv >> `basename $0`
          rm -f headerless_report.tsv
          head `basename $0`
    out:
    - output_file

  feature_select:
    run: ../tools/feature-select-sql.cwl
    in:
      feature_file: rename_header/output_file
      sql_query: sql_query
      columns:
        source: columns
        valueFrom: $("DISTINCT " + self.join(", "))   # multiple peaks can have the same coordinates but different abssummit, so we need to use DISTINCT
      header: header
    out:
    - filtered_file
    - stdout_log
    - stderr_log

  add_header:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
      - class: ScatterFeatureRequirement
      - class: ShellCommandRequirement
      inputs:
        script:
          type: string?
          default: |
            printf "Chrom\tStart\tEnd\tName\tScore\tStrand\n" > genelist-filtered-set-w-header.bed
            cat $0 >> genelist-filtered-set-w-header.bed
          inputBinding:
            position: 1
        headerless_bed:
          type: File
          inputBinding:
            position: 2
      outputs:
        filtered_file_with_header:
          type: File
          outputBinding:
            glob: genelist-filtered-set-w-header.bed
      baseCommand: ["bash", "-c"]
    in:
      headerless_bed: feature_select/filtered_file
    out:
    - filtered_file_with_header


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Filter differentially bound sites for heatmap analysis"
label: "Filter differentially bound sites for heatmap analysis"
s:alternateName: "Filter differentially bound sites from DiffBind analysis to be used with deepTools heatmap analysis"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/filter-diffbind-for-heatmap.cwl
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
  Filter DiffBind results for deepTools heatmap analysis
  ======================================================

  Filter differentially bound sites from DiffBind analysis to be used with deepTools heatmap analysis