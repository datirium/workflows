cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/deseq:v0.0.7
inputs:
  expression_files:
    type: File[]
    inputBinding:
      prefix: --expression
    doc: |
      Path to the TSV/CSV files with expression data.
      All files should have the following header:
      RefseqId GeneId Chrom TxStart TxEnd Strand TotalReads Rpkm
  expression_names:
    type: string[]
    inputBinding:
      prefix: --aliases
    doc: |
      Unique names for files provided in --expression,
      no special characters or spaces are allowed.
      Number and order of the names should corresponds
      to values from --expression.
  metadata_file:
    type: File
    inputBinding:
      prefix: --metadata
    doc: |
      Path to the TSV/CSV file to provide metadata for the
      samples from --expression. First column should have
      the name 'sample', other columns may have arbitrary names.
      The values from the 'sample' column should correspond to
      the values provided in --aliases. For a proper --contrast
      intepretation, values defined in each column should not be
      used in other columns. All metadata columns are treated as
      factors (no covariates are supported).
  design_formula:
    type: string
    inputBinding:
      prefix: --design
    doc: |
      Design formula. Should start with ~ and include terms from
      the --metadata table.
  reduced_formula:
    type: string?
    inputBinding:
      prefix: --reduced
    doc: |
      Reduced formula with the term(s) of interest removed.
      Should start with ~. If provided, force DESeq2 to run
      LRT test instead of the Wald.
  contrast:
    type: string?
    inputBinding:
      prefix: --contrast
    doc: |
      Contrast to be be applied for the output, formatted as
      a mathematical formula of values from the --metadata table.
      If not provided, all possible combinations of values from
      the metadata columns present in the --design will be used
      (results will be merged giving the priority to significantly
      differentially expressed genes with higher absolute
      log2FoldChange values).
  base:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      prefix: --base
    doc: |
      Value(s) from each metadata file column(s) to be set as
      the base level(s). Number and order of provided values should
      correspond the order of columns in --metadata file. Default:
      define base levels alphabetically for each metadata column.
  feature_type:
    type:
    - 'null'
    - type: enum
      symbols:
      - gene
      - transcript
    inputBinding:
      prefix: --type
    doc: |
      Feature type to use for differential expression.
      If set to 'gene', use 'GeneId' column from the provided in --expression files.
      If set to 'transcript', use 'RefseqId' from the provided in --expression files.
      Default: gene
  excluded_features:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      prefix: --exclude
    doc: |
      Features to be excluded from the differential expression analysis.
      Default: include all features
  normalization_method:
    type:
    - 'null'
    - type: enum
      symbols:
      - vst
      - rlog
    inputBinding:
      prefix: --norm
    doc: |
      Read counts normalization for the exploratory visualization analysis.
      Use 'vst' for medium-to-large datasets (n > 30) and 'rlog' for
      small datasets (n < 30), when there is a wide range of sequencing
      depth across samples.
      Default: vst
  remove:
    type: string?
    inputBinding:
      prefix: --remove
    doc: |
      Column from the metadata file to remove batch effect
      before running differential expression analysis. If
      present, all components that include this term will be
      removed from the design and reduced formulas.
      Default: do not remove batch effect
  cluster_method:
    type:
    - 'null'
    - type: enum
      symbols:
      - row
      - column
      - both
    inputBinding:
      prefix: --cluster
    doc: |
      Hopach clustering method to be run on normalized read counts for the
      exploratory visualization analysis. Default: do not run clustering
  row_distance:
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
    doc: |
      Distance metric for HOPACH row clustering. Ignored if --cluster is not
      provided. Default: cosangle
  column_distance:
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
      prefix: --columndist
    doc: |
      Distance metric for HOPACH column clustering. Ignored if --cluster is not
      provided. Default: euclid
  center_row:
    type: boolean?
    inputBinding:
      prefix: --center
    doc: |
      Apply mean centering for feature expression prior to running
      clustering by row. Ignored when --cluster is not row or both.
      Default: do not centered
  selected_features:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      prefix: --label
    doc: |
      Features of interest to label on the generated volcanot plot. Default:
      top 10 features with the highest and the lowest log2 fold change
      expression values.
  maximum_padj:
    type: float?
    inputBinding:
      prefix: --padj
    doc: |
      In the exploratory visualization analysis output only features with
      adjusted P-value not bigger than this value. Default: 0.05
  minimum_logfc:
    type: float?
    inputBinding:
      prefix: --logfc
    doc: |
      In the exploratory visualization analysis output only features with
      absolute log2FoldChange bigger or equal to this value. Default: 0
  export_pdf_plots:
    type: boolean?
    inputBinding:
      prefix: --pdf
    doc: |
      Export plots in PDF.
      Default: false
  output_prefix:
    type: string?
    inputBinding:
      prefix: --output
    doc: |
      Output prefix for generated files
  threads:
    type: int?
    inputBinding:
      prefix: --cpus
    doc: |
      Number of cores/cpus to use. Default: 1
outputs:
  diff_expr_features:
    type: File
    outputBinding:
      glob: '*_diff_expr_features.tsv'
    doc: |
      TSV file with not filtered differentially expressed features
  read_counts_gct:
    type: File
    outputBinding:
      glob: '*_norm_read_counts.gct'
    doc: |
      GCT file with normalized, optionally batch corrected, read counts
  volcano_plot_png:
    type: File?
    outputBinding:
      glob: '*_volcano_plot.png'
    doc: |
      Volcano plot of differentially expressed features.
      PNG format
  volcano_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_volcano_plot.pdf'
    doc: |
      Volcano plot of differentially expressed features.
      PDF format
  pca_plot_png:
    type: File?
    outputBinding:
      glob: '*_pca_plot.png'
    doc: |
      PCA plot of normalized, optionally batch corrected,
      read counts based on the top 500 features selected
      by the highest row variance
      PNG format
  pca_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_pca_plot.pdf'
    doc: |
      PCA plot of normalized, optionally batch corrected,
      read counts based on the top 500 features selected
      by the highest row variance
      PDF format
  mds_plot_html:
    type: File?
    outputBinding:
      glob: '*_mds_plot.html'
    doc: |
      MDS plot of normalized, optionally batch corrected,
      read counts
      HTML format
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- run_deseq_manual.R
stdout: deseq_manual_factor_stdout.log
stderr: deseq_manual_factor_stderr.log
label: DESeq2 Multi-factor Analysis
doc: |
  DESeq2 Multi-factor Analysis

  Runs DeSeq2 multi-factor analysis with manual control over major parameters
