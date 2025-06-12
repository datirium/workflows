cwlVersion: v1.0
class: CommandLineTool
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/hopach:v0.0.10
inputs:
  expression_files:
    type: File[]
    inputBinding:
      prefix: --input
    doc: Input CSV/TSV files with RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand, TotalReads, Rpkm columns
  expression_aliases:
    type:
    - 'null'
    - string[]
    inputBinding:
      prefix: --name
    doc: 'Input aliases, the order corresponds to --input order. Default: basename of --input files'
  genelist_file:
    type: File?
    inputBinding:
      prefix: --genelist
    doc: Filter genes by the list from the file. Headerless, 1 gene per line
  target_column:
    type: string?
    inputBinding:
      prefix: --target
    doc: 'Target column to be used by hopach clustering. Default: Rpkm'
  combine:
    type:
    - 'null'
    - string[]
    inputBinding:
      prefix: --combine
    doc: 'Combine inputs by columns names. Default: RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand'
  cluster_method:
    type:
    - 'null'
    - type: enum
      symbols:
      - row
      - column
      - both
    inputBinding:
      prefix: --method
    doc: 'Cluster method. Default: both'
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
    inputBinding:
      prefix: --rowdist
    doc: 'Distance metric for row clustering. Default: cosangle'
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
    inputBinding:
      prefix: --coldist
    doc: 'Distance metric for column clustering. Default: euclid'
  row_logtransform:
    type: boolean?
    inputBinding:
      prefix: --rowlogtransform
    doc: 'Log2 transform input data prior to running row clustering. Default: false'
  col_logtransform:
    type: boolean?
    inputBinding:
      prefix: --collogtransform
    doc: 'Log2 transform input data prior to running column clustering. Default: false'
  row_center:
    type:
    - 'null'
    - type: enum
      symbols:
      - mean
      - median
    inputBinding:
      prefix: --rowcenter
    doc: 'Center rows prior to running row clustering. Default: not centered'
  col_center:
    type:
    - 'null'
    - type: enum
      symbols:
      - mean
      - median
    inputBinding:
      prefix: --colcenter
    doc: 'Center columns prior to running column clustering. Default: not centered'
  row_normalize:
    type: boolean?
    inputBinding:
      prefix: --rownorm
    doc: 'Normalize rows prior to running row clustering. Default: not normalized'
  col_normalize:
    type: boolean?
    inputBinding:
      prefix: --colnorm
    doc: 'Normalize columns prior to running column clustering. Default: not normalized'
  row_min:
    type: float?
    inputBinding:
      prefix: --rowmin
    doc: 'Exclude rows from clustering by the min value of a target column. Default: 0'
  keep_discarded:
    type: boolean?
    inputBinding:
      prefix: --rowkeep
    doc: 'Append excluded rows to the output table after clustering is finished. Default: false'
  palette:
    type:
    - 'null'
    - string[]
    inputBinding:
      prefix: --palette
    doc: 'Palette color names. Default: red, black, green'
  output_prefix:
    type: string?
    inputBinding:
      prefix: --output
    doc: 'Output prefix. Default: hopach'
outputs:
  clustering_results:
    type: File
    outputBinding:
      glob: '*_clustering.tsv'
    doc: Hopach clustering results
  clustering_matrix:
    type: File
    outputBinding:
      glob: '*_tree.cdt'
    doc: Hopach clustering results as cdt file
  clustering_row_dendogram:
    type: File?
    outputBinding:
      glob: '*_tree.gtr'
    doc: Hopach row clustering dendogram as gtr file
  clustering_column_dendogram:
    type: File?
    outputBinding:
      glob: '*_tree.atr'
    doc: Hopach column clustering dendogram as atr file
  heatmap_png:
    type: File
    outputBinding:
      glob: '*_heatmap.png'
    doc: Heatmap ordered by hopach clustering results
  column_clustering_labels:
    type: File?
    outputBinding:
      glob: '*_column_clustering_labels.tsv'
    doc: Hopach column clustering labels
  row_distance_matrix_png:
    type: File?
    outputBinding:
      glob: '*_row_dist_matrix.png'
    doc: Row distance matrix
  col_distance_matrix_png:
    type: File?
    outputBinding:
      glob: '*_column_dist_matrix.png'
    doc: Column distance matrix
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- hopach_order.R
stderr: hopach_stderr.log
stdout: hopach_stdout.log
label: HOPACH - Hierarchical Ordered Partitioning and Collapsing Hybrid
doc: |
  Runs hopach clustering algorithm with the combined by specific columns input files.
  Works with minimum two genelist files.
