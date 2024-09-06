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
    - "filter-peaks-for-heatmap.cwl"
    - "genelists-sets.cwl"
  samples_rnaseq:
    - "mirna-mirdeep2-se.cwl"
    - "rnaseq-se.cwl"
    - "rnaseq-pe.cwl"
    - "rnaseq-se-dutp.cwl"
    - "rnaseq-pe-dutp.cwl"
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-pe-smarter-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"
    - "trim-quantseq-mrnaseq-se-strand-specific.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  genelist_names:
    type:
    - "null"
    - string[]
    default: null
    label: "Genelist sample names"
    doc: "Array of genelist aliases/sample names, param (-a)"
    'sd:upstreamSource': "genelists/alias"
    'sd:localLabel': true

  genelist_feature_files:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3475"
    label: "feature files from genelist inputs"
    doc: "Array of TSV files with differential genes from DESeq or diffbind pipelines, param (-b)"
    'sd:upstreamSource': "genelists/feature_file"

  genelist_filtered_files:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3475"
    label: "Filtered differential genes"
    doc: "Array of filtered differential genelists from DESeq or diffbind pipelines, param (-c)"
    'sd:upstreamSource': "genelists/filtered_file"

  sample_names_rnaseq:
    type:
      - "null"
      - string[]?
    default: null
    label: "Sample names for RNA-Seq experiments"
    doc: "Array of aliases for RNA-Seq experiments for column metadata, param (-e)"
    'sd:upstreamSource': "samples_rnaseq/alias"
    'sd:localLabel': true

  datafiles_rnaseq:
    type:
    - "null"
    - File[]?
    default: null
    format: "http://edamontology.org/format_3475"
    label: "Array of sample TSV files containing gene annotations with associated TotalReads and Rpkm counts"
    doc: "Array of sample TSV files containing gene annotations with associated TotalReads and Rpkm counts, param (-g)"
    'sd:upstreamSource': "samples_rnaseq/rpkm_genes"

  threads:
    type: int?
    default: 1
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading."
    'sd:layout':
      advanced: true


outputs:

  master_samplesheet_scaled:
    type: File
    label: "contains formatted information of the input data and files for scaled heatmaps"
    format: "http://edamontology.org/format_3475"
    doc: "contains formatted information of the input data and files for scaled heatmaps"
    outputSource: data_integration/master_samplesheet_scaled

  master_samplesheet_vst:
    type: File
    label: "contains formatted information of the input data and files for vst heatmap"
    format: "http://edamontology.org/format_3475"
    doc: "contains formatted information of the input data and files for vst heatmap"
    outputSource: data_integration/master_samplesheet_vst

  master_samplesheet_vst_zscore:
    type: File
    label: "contains formatted information of the input data and files for vst zscore heatmap"
    format: "http://edamontology.org/format_3475"
    doc: "contains formatted information of the input data and files for vst zscore heatmap"
    outputSource: data_integration/master_samplesheet_vst_zscore

  heatmap_TotalReads_html:
    type: File
    outputSource: data_integration/heatmap_TotalReads_html
    label: "Heatmap TotalReads (raw total read counts)"
    doc: |
      html of morpheus heatmap with preconfigured settings, TotalReads, no data scaling
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  heatmap_vst_html:
    type: File
    outputSource: data_integration/heatmap_vst_html
    label: "Heatmap VST (VST normalized TotalReads)"
    doc: |
      html of morpheus heatmap with preconfigured settings, VST values, no data scaling
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  heatmap_vst_zscore_html:
    type: File
    outputSource: data_integration/heatmap_vst_zscore_html
    label: "Heatmap VST Z-scores (Z-score of VST)"
    doc: |
      html of morpheus heatmap with preconfigured settings, Z-scores of VST values, no data scaling
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  heatmap_Rpkm_html:
    type: File
    outputSource: data_integration/heatmap_Rpkm_html
    label: "Heatmap RPKM (RPKM normalized TotalReads)"
    doc: |
      html of morpheus heatmap with preconfigured settings, RPKM, no data scaling
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  heatmap_scaled100_html:
    type: File
    outputSource: data_integration/heatmap_scaled100_html
    label: "Heatmap RPKM (scaled 0-99, no percentile cutoff)"
    doc: |
      html of morpheus heatmap with preconfigured settings, RPKM, data scaled 0-99, no percentile cutoff
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  heatmap_scaled99_html:
    type: File
    outputSource: data_integration/heatmap_scaled99_html
    label: "Heatmap RPKM (scaled 0-99, max value set to 99th percentile RPKM value)"
    doc: |
      html of morpheus heatmap with preconfigured settings, RPKM, data scaled 0-99, max value set to 99th percentile RPKM value
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  heatmap_scaled95_html:
    type: File
    outputSource: data_integration/heatmap_scaled95_html
    label: "Heatmap RPKM (scaled 0-99, max value set to 95th percentile RPKM value)"
    doc: |
      html of morpheus heatmap with preconfigured settings, RPKM, data scaled  0-99, max value set to 95th percentile RPKM value
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  log_file_stdout:
    type: File
    outputSource: data_integration/log_file_stdout
    label: "log file stdout"

  log_file_stderr:
    type: File
    outputSource: data_integration/log_file_stderr
    label: "log file stderr"


steps:

  data_integration:
    run: ../tools/genelists-deseq-only.cwl
    in:
      threads: threads
      genelist_names: genelist_names
      feature_files: genelist_feature_files
      filtered_files: genelist_filtered_files
      sample_names_rnaseq: sample_names_rnaseq
      expression_files: datafiles_rnaseq
    out:
      - master_samplesheet_scaled
      - master_samplesheet_vst
      - master_samplesheet_vst_zscore
      - heatmap_TotalReads_html
      - heatmap_vst_html
      - heatmap_vst_zscore_html
      - heatmap_Rpkm_html
      - heatmap_scaled100_html
      - heatmap_scaled99_html
      - heatmap_scaled95_html
      - log_file_stdout
      - log_file_stderr


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Genelists heatmap - RNA-seq expression data visualized"
label: "Genelists heatmap - RNA-seq expression data visualized"
s:alternateName: "Genelists heatmap with RNA-Seq expression data only"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/genelists-deseq-only.cwl
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
  # Genelists heatmap - RNA-seq expression data visualized

  This visualization workflow takes as input 1 or more genelists derived from the DESeq and/or diffbind workflows along with user-selected samples and visualizes RNA-Seq expression data in a single morpheus heatmap.

  ### __References__
    - Morpheus, https://software.broadinstitute.org/morpheus