cwlVersion: v1.0
class: Workflow
requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
- class: MultipleInputFeatureRequirement
sd:upstream:
  rnaseq_sample:
  - rnaseq-se.cwl
  - rnaseq-pe.cwl
  - rnaseq-se-dutp.cwl
  - rnaseq-pe-dutp.cwl
  - rnaseq-se-dutp-mitochondrial.cwl
  - rnaseq-pe-dutp-mitochondrial.cwl
  - trim-rnaseq-pe.cwl
  - trim-rnaseq-se.cwl
  - trim-rnaseq-pe-dutp.cwl
  - trim-rnaseq-se-dutp.cwl
inputs:
  alias:
    type: string
    label: Experiment short name/Alias
    sd:preview:
      position: 1
  expression_files:
    type: File[]
    format: http://edamontology.org/format_3752
    label: Isoform expression files
    doc: Isoform expression files
    sd:upstreamSource: rnaseq_sample/rpkm_isoforms
    sd:localLabel: true
  expression_aliases:
    type:
    - 'null'
    - string[]
    label: Isoform expression file aliases
    doc: Aliases to make the legend for generated plots. Order corresponds to the isoform expression files
    sd:upstreamSource: rnaseq_sample/alias
  genelist_file:
    type: File?
    format: http://edamontology.org/format_2330
    label: Gene list to filter clustered genes. Headerless TSV/CSV file with 1 gene per line
    doc: Gene list to filter clustered genes. Headerless TSV/CSV file with 1 gene per line
  group_by:
    type:
    - 'null'
    - type: enum
      symbols:
      - isoforms
      - genes
      - common tss
    default: genes
    label: Group by
    doc: 'Grouping method for features: isoforms, genes or common tss'
  target_column:
    type:
    - 'null'
    - type: enum
      symbols:
      - Rpkm
      - TotalReads
    default: Rpkm
    label: Target column
    doc: Target column name from expression files to be used by hopach
  row_min:
    type: float?
    default: 0
    label: Min value for target column
    doc: Min value for target column
  keep_discarded:
    type: boolean?
    default: false
    label: Keep discarded rows
    doc: Keep discarded by threshold parameter rows at the end of the output file
  cluster_method:
    type:
    - 'null'
    - type: enum
      symbols:
      - row
      - column
      - both
    default: both
    label: Cluster method
    doc: Cluster method
  row_dist_metric:
    type:
    - 'null'
    - type: enum
      symbols:
      - cosangle
      - abscosangle
      - euclid
      - abseuclid
      - cor
      - abscor
    default: cosangle
    label: Distance metric for row clustering
    doc: Algorithm to be used for distance matrix calculation before running hopach row clustering
    sd:layout:
      advanced: true
  row_logtransform:
    type: boolean?
    default: false
    label: Log2 row transform
    doc: Log2 transform input data to prior running hopach row clustering
    sd:layout:
      advanced: true
  row_center:
    type:
    - 'null'
    - type: enum
      symbols:
      - mean
      - median
    default: mean
    label: Center row values
    doc: Center row values
    sd:layout:
      advanced: true
  row_normalize:
    type: boolean?
    default: true
    label: Normalize row values
    doc: Normalize row values
    sd:layout:
      advanced: true
  col_dist_metric:
    type:
    - 'null'
    - type: enum
      symbols:
      - cosangle
      - abscosangle
      - euclid
      - abseuclid
      - cor
      - abscor
    default: euclid
    label: Distance metric for column clustering
    doc: Algorithm to be used for distance matrix calculation before running hopach column clustering
    sd:layout:
      advanced: true
  col_logtransform:
    type: boolean?
    default: false
    label: Log2 column transform
    doc: Log2 transform input data to prior running hopach column clustering
    sd:layout:
      advanced: true
  col_center:
    type:
    - 'null'
    - type: enum
      symbols:
      - mean
      - median
    default: mean
    label: Center column values
    doc: Center column values
    sd:layout:
      advanced: true
  col_normalize:
    type: boolean?
    default: true
    label: Normalize column values
    doc: Normalize column values
    sd:layout:
      advanced: true
  palette:
    type:
    - 'null'
    - string[]
    default:
    - black
    - red
    - yellow
    label: Custom palette color list
    doc: Color list for custom palette
    sd:layout:
      advanced: true
outputs:
  clustering_results:
    type: File
    format: http://edamontology.org/format_3475
    label: Combined clustered expression file
    doc: Combined by RefseqId, GeneId, Chrom, TxStart, TxEnd and Strand clustered expression file
    outputSource: hopach/clustering_results
    sd:visualPlugins:
    - syncfusiongrid:
        tab: Hopach Clustering Results
        Title: Combined clustered expression file
  clustering_matrix:
    type: File
    format: http://edamontology.org/format_3475
    label: Hopach clustering results as cdt file
    doc: Hopach clustering results as cdt file
    outputSource: hopach/clustering_matrix
  clustering_row_dendogram:
    type: File?
    format: http://edamontology.org/format_3475
    label: Hopach row clustering dendogram as gtr file
    doc: Hopach row clustering dendogram as gtr file
    outputSource: hopach/clustering_row_dendogram
  clustering_column_dendogram:
    type: File?
    format: http://edamontology.org/format_3475
    label: Hopach column clustering dendogram as atr file
    doc: Hopach column clustering dendogram as atr file
    outputSource: hopach/clustering_column_dendogram
  column_clustering_labels:
    type: File?
    format: http://edamontology.org/format_3475
    label: Column cluster labels
    doc: Column cluster labels
    outputSource: hopach/column_clustering_labels
    sd:visualPlugins:
    - syncfusiongrid:
        tab: Hopach Clustering Results
        Title: Column cluster labels
  heatmap_png:
    type: File
    label: Heatmap
    format: http://edamontology.org/format_3603
    doc: Heatmap plot
    outputSource: hopach/heatmap_png
    sd:visualPlugins:
    - image:
        tab: Plots
        Caption: Heatmap
  row_distance_matrix_png:
    type: File?
    label: Row Distance Matrix
    format: http://edamontology.org/format_3603
    doc: Row distance matrix plot. Clusters of similar features will appear as blocks on the diagonal of the matrix
    outputSource: hopach/row_distance_matrix_png
    sd:visualPlugins:
    - image:
        tab: Plots
        Caption: Row Distance Matrix
  col_distance_matrix_png:
    type: File?
    label: Column Distance Matrix
    format: http://edamontology.org/format_3603
    doc: Column distance matrix plot. Clusters of similar features will appear as blocks on the diagonal of the matrix
    outputSource: hopach/col_distance_matrix_png
    sd:visualPlugins:
    - image:
        tab: Plots
        Caption: Column Distance Matrix
  hopach_stdout_log:
    type: File
    format: http://edamontology.org/format_2330
    label: HOPACH stdout log
    doc: HOPACH stdout log
    outputSource: hopach/stdout_log
  hopach_stderr_log:
    type: File
    format: http://edamontology.org/format_2330
    label: HOPACH stderr log
    doc: HOPACH stderr log
    outputSource: hopach/stderr_log
steps:
  group_isoforms:
    run: ../tools/group-isoforms-batch.cwl
    in:
      isoforms_file: expression_files
    out:
    - genes_file
    - common_tss_file
  hopach:
    run: ../tools/hopach.cwl
    in:
      expression_files:
        source:
        - group_by
        - expression_files
        - group_isoforms/genes_file
        - group_isoforms/common_tss_file
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
      expression_aliases: expression_aliases
      genelist_file: genelist_file
      target_column: target_column
      cluster_method: cluster_method
      row_dist_metric: row_dist_metric
      col_dist_metric: col_dist_metric
      row_logtransform: row_logtransform
      col_logtransform: col_logtransform
      row_center: row_center
      col_center: col_center
      row_normalize: row_normalize
      col_normalize: col_normalize
      row_min: row_min
      keep_discarded: keep_discarded
      palette: palette
    out:
    - clustering_results
    - clustering_matrix
    - clustering_row_dendogram
    - clustering_column_dendogram
    - column_clustering_labels
    - heatmap_png
    - row_distance_matrix_png
    - col_distance_matrix_png
    - stdout_log
    - stderr_log
label: HOPACH - Hierarchical Ordered Partitioning and Collapsing Hybrid
doc: |
  Hierarchical Ordered Partitioning and Collapsing Hybrid (HOPACH)
  ===============================================================

  The HOPACH clustering algorithm builds a hierarchical tree of clusters by recursively partitioning a data set, while
  ordering and possibly collapsing clusters at each level. The algorithm uses the Mean/Median Split Silhouette (MSS) criteria
  to identify the level of the tree with maximally homogeneous clusters. It also runs the tree down to produce a final
  ordered list of the elements. The non-parametric bootstrap allows one to estimate the probability that each element
  belongs to each cluster (fuzzy clustering).
sd:version: 100
