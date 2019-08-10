cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  rnaseq_sample:
    - "rnaseq-se.cwl"
    - "rnaseq-pe.cwl"
    - "rnaseq-se-dutp.cwl"
    - "rnaseq-pe-dutp.cwl"
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  expression_files:
    type: File[]
    format: "http://edamontology.org/format_3752"
    label: "Isoform expression files"
    doc: "Isoform expression files"
    'sd:upstreamSource': "rnaseq_sample/rpkm_isoforms"
    'sd:localLabel': true

  legend_names:
    type:
      - "null"
      - string[]
    label: "Isoform expression file aliases"
    doc: "Aliases to make the legend for generated plots. Order corresponds to the isoform expression files"
    # 'sd:upstreamSource': "rnaseq_sample/alias"
    # 'sd:localLabel': true

  group_by:
    type:
      - "null"
      - type: enum
        symbols: ["isoforms", "genes", "common tss"]
    default: "genes"
    label: "Group by"
    doc: "Grouping method for features: isoforms, genes or common tss"

  target_column:
    type:
      - "null"
      - type: enum
        symbols: ["Rpkm", "TotalReads"]
    default: "Rpkm"
    label: "Target column"
    doc: "Target column name from expression files to be used by Hopach"

  threshold:
    type: float?
    default: 0
    label: "Target column threshold value"
    doc: "Min value for target column"

  keep_discarded:
    type: boolean?
    default: false
    label: "Append discarded by threshold rows at the bottom"
    doc: "Keep discarded by threshold parameter rows at the end of the output file"

  dist_metric:
    type:
      - "null"
      - type: enum
        symbols: ["cosangle", "abscosangle", "euclid", "abseuclid", "cor", "abscor"]
    default: "cosangle"
    label: "Distance metric"
    doc: "Algorithm to be used for distance matrix calculation before running hopach clustering"

  logtransform:
    type: boolean?
    default: false
    label: "Log2 transform"
    doc: "Log2 transform input data prior running hopach clustering"

  export_heatmap:
    type: boolean?
    default: false
    label: "Export heatmap"
    doc: "Export heatmap plot to png"
    'sd:layout':
      advanced: true

  export_distance_matrix:
    type: boolean?
    default: false
    label: "Export distance matrix plot"
    doc: "Export distance matrix plot to png"
    'sd:layout':
      advanced: true

  export_variability_plot:
    type: boolean?
    default: false
    label: "Export cluster variability plot"
    doc: "Export cluster variability plot"
    'sd:layout':
      advanced: true

  palette:
    type:
      - "null"
      - string[]
    default: ["black", "red", "yellow"]
    label: "Custom palette color list"
    doc: "Color list for custom palette"
    'sd:layout':
      advanced: true


outputs:

  ordered_genelist:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Combined expression file ordered by hopach clustering results"
    doc: "Combined by RefseqId, GeneId, Chrom, TxStart, TxEnd and Strand expression file. Rows order correspond to the Hopach clustering results"
    outputSource: hopach/ordered_genelist
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Hopach Clustering Analysis'
        Title: 'Hopach ordered expression data'

  distance_matrix_png:
    type: File
    label: "Distance Matrix"
    format: "http://edamontology.org/format_3603"
    doc: "Distance matrix plot. Clusters of similar features will appear as blocks on the diagonal of the matrix"
    outputSource: hopach/distance_matrix_png
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'Distance Matrix'

  heatmap_png:
    type: File
    label: "Heatmap"
    format: "http://edamontology.org/format_3603"
    doc: "Heatmap plot. Row ordering corresponds to the Hopach clustering results"
    outputSource: hopach/heatmap_png
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'Heatmap'

  variability_plot_png:
    type: File
    label: "Cluster Variability"
    format: "http://edamontology.org/format_3603"
    doc: |
      Cluster variability plot. Every horizontal bar represents a feature.
      If the bar is all or mostly one color, then the feature is estimated to
      belong strongly to that cluster
    outputSource: hopach/variability_plot_png
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'Cluster Variability'

steps:

  group_isoforms:
    run: ../subworkflows/group-isoforms-batch.cwl
    in:
      isoforms_file: expression_files
    out:
      - genes_file
      - common_tss_file

  hopach:
    run: ../tools/hopach.cwl
    in:
      genelist_files:
        source: [group_by, expression_files, group_isoforms/genes_file, group_isoforms/common_tss_file]
        valueFrom: |
          ${
              if (self[0] == "isoforms") {
                return self[1];
              } else if (self[0] == "genes") {
                return self[2];
              } else {
                return self[3];
              }
          }
      export_heatmap: export_heatmap
      export_distance_matrix: export_distance_matrix
      export_variability_plot: export_variability_plot
      legend_names: legend_names
      target_column: target_column
      dist_metric: dist_metric
      logtransform: logtransform
      keep_discarded: keep_discarded
      threshold: threshold
      palette: palette
    out:
      - ordered_genelist
      - distance_matrix_png
      - heatmap_png
      - variability_plot_png


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "HOPACH - Hierarchical Ordered Partitioning and Collapsing Hybrid"
label: "HOPACH - Hierarchical Ordered Partitioning and Collapsing Hybrid"
s:alternateName: "The HOPACH clustering algorithm builds a hierarchical tree of clusters by recursively partitioning a data set"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/hopach.cwl
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
          s:email: mailto:michael.kotliar@cchmc.org
          s:sameAs:
          - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Hierarchical Ordered Partitioning and Collapsing Hybrid (HOPACH)
  ===============================================================

  The HOPACH clustering algorithm builds a hierarchical tree of clusters by recursively partitioning a data set, while
  ordering and possibly collapsing clusters at each level. The algorithm uses the Mean/Median Split Silhouette (MSS) criteria
  to identify the level of the tree with maximally homogeneous clusters. It also runs the tree down to produce a final
  ordered list of the elements. The non-parametric bootstrap allows one to estimate the probability that each element
  belongs to each cluster (fuzzy clustering).
