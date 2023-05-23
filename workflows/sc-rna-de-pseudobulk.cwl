cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var split_features = function(line) {
          function get_unique(value, index, self) {
            return self.indexOf(value) === index && value != "";
          }
          var splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };


'sd:upstream':
  sc_tools_sample:
  - "sc-rna-cluster.cwl"
  - "sc-ctype-assign.cwl"
  - "sc-wnn-cluster.cwl"
  - "sc-rna-da-cells.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Single-cell Cluster or Manual Cell Type Assignment Analysis"
    doc: |
      Single-cell analysis run through the
      clustering or cell type assignment
      pipelines.
    'sd:upstreamSource': "sc_tools_sample/seurat_data_rds"
    'sd:localLabel': true

  datasets_metadata:
    type: File?
    label: "TSV/CSV file to assign categories per sample"
    doc: |
      If selected single-cell analysis was run
      with the data aggregated from multiple
      samples, you can optionally provide tab-
      delimited or comma-separated file to
      assign additional categories per sample.
      First column should be named 'library_id'
      and include all sample names from the
      selected single-cell analysis regardless
      whether filtering by barcodes was applied
      or not. All other columns may have
      arbitrary names.

  barcodes_data:
    type: File?
    label: "TSV/CSV file to filter cells by barcodes"
    doc: |
      Loaded single-cell data can be optionally
      prefiltered by selected cell barcodes.
      Provided tab-delimited or comma-separated
      file should have the first column named
      'barcode'. If this file includes any other
      columns, they will be used to assign
      additional categories per cell.

  groupby:
    type: string?
    default: null
    label: "Category to group cells for optional subsetting"
    doc: |
      Before running differential expression
      analysis input data can be optionally
      prefiltered to include only certain
      values from the specific category.
      Here we define the name of that
      category.

  subset:
    type: string?
    default: null
    label: "List of values to subset cells from the selected category"
    doc: |
      If the category to group cells for
      optional subsetting was provided,
      here we define which values should
      be included into analysis.

  splitby:
    type: string
    label: "Category to split cell into two groups"
    doc: |
      All remaining after optional prefiltering
      steps cells will be split into two groups
      for gene expression comparison.

  first_cond:
    type: string
    label: "Value from the selected category to define the first group of cells"
    doc: |
      Cells for which the selected category
      includes provided value will be used
      as the first group for differential
      expression comparison. Direction of
      comparison is second vs first groups.

  second_cond:
    type: string
    label: "Value from the selected category to define the second group of cells"
    doc: |
      Cells for which the selected category
      includes provided value will be used
      as the second group for differential
      expression comparison. Direction of
      comparison is second vs first groups.

  analysis_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "wilcoxon (by cells, no batches)"                   # (wilcox) Wilcoxon Rank Sum test
      - "likelihood-ratio (by cells, no batches)"           # (bimod) Likelihood-ratio test
      - "t-test (by cells, no batches)"                     # (t) Student's t-test
      - "negative-binomial (by cells, models batches)"      # (negbinom) Negative Binomial Generalized Linear Model (supports --batchby)
      - "poisson (by cells, models batches)"                # (poisson) Poisson Generalized Linear Model (supports --batchby)
      - "logistic-regression (by cells, models batches)"    # (LR) Logistic Regression (supports --batchby)
      - "mast (by cells, models batches)"                   # (MAST) MAST package (supports --batchby)
      - "deseq (pseudo bulk, models batches)"               # DESeq2 Wald test on pseudobulk aggregated gene expression
      - "deseq-lrt (pseudo bulk, models batches)"           # DESeq2 LRT test on pseudobulk aggregated gene expression
    default: wilcoxon
    label: "Test type to use in differential expression analysis"
    doc: |
      Test type to use in the differential
      expression analysis. If set to deseq
      or deseq-lrt, gene expression will be
      aggregated to the pseudobulk form per
      sample. Othwerwise, analysis will be
      run on the cells level. If deseq is
      selected, the pair-wise Wald test will
      be used. For deseq-lrt, the Likelihood
      Ratio Test will be applied between
      design and reduced formulas. The
      reduced formula will look like ~1 if
      grouping by batches is omitted or will
      be set to the category defined as
      batches.

  batchby:
    type: string?
    default: null
    label: "Category to model batch effect"
    doc: |
      If selected test type supports batch
      effect modeling, the provided category
      will be used to group cells into
      batches. For deseq and deseq-lrt tests
      batch modeling will result in adding it
      into the design formula. For negative-
      binomial, poisson, logistic-regression,
      or mast tests grouping by batches will
      be used as a latent variable in the
      FindMarkers function.

  maximum_padj:
    type: float?
    default: 0.05
    label: "Maximum adjusted P-value for genes displayed on the heatmap"
    doc: |
      When generating gene expression heatmap
      per cell output only differentially
      expressed genes with the adjusted P-value
      not bigger than this value.

  genes_of_interest:
    type: string?
    default: null
    label: "Genes of interest to be shown on the plots"
    doc: |
      Genes of interest to be shown on the
      volcano, violin, and UMAP plots.
    'sd:layout':
      advanced: true

  exclude_pattern:
    type: string?
    default: null
    label: "Regex pattern to identify and exclude specific genes from the analysis"
    doc: |
      Regex pattern to identify and exclude
      specific genes from the differential
      expression analysis (not case-sensitive).
      If any of such genes were selected as
      genes of interest to be shown on the plots,
      they will be excluded from there as well.
    'sd:layout':
      advanced: true

  cluster_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "row"
      - "column"
      - "both"
      - "none"
    default: "row"
    label: "Clustering method for gene expression data"
    doc: |
      Clustering method to be run on
      the normalized read counts data.
      "column" and "both" options are
      supported only when using deseq
      or desey-lrt tests for which gene
      expression data aggregated to the
      pseudobulk form.
    'sd:layout':
      advanced: true

  row_distance:
    type:
    - "null"
    - type: enum
      symbols:
      - "cosangle"
      - "abscosangle"
      - "euclid"
      - "abseuclid"
      - "cor"
      - "abscor"
    default: "cosangle"
    label: "Distance metric for row clustering"
    doc: |
      Distance metric for row clustering.
      Ignored if clustering method is set
      to "column" or "none".
    'sd:layout':
      advanced: true

  column_distance:
    type:
    - "null"
    - type: enum
      symbols:
      - "cosangle"
      - "abscosangle"
      - "euclid"
      - "abseuclid"
      - "cor"
      - "abscor"
    default: "euclid"
    label: "Distance metric for column clustering"
    doc: |
      Distance metric for column clustering.
      Ignored if clustering method is set
      to "row" or "none".
    'sd:layout':
      advanced: true

  center_row:
    type: boolean?
    default: true
    label: "Gene expression mean centering for clustering by row"
    doc: |
      Apply mean centering for gene
      expression prior to running
      clustering by row. Ignored if
      clustering method is set to
      "column" or "none".
    'sd:layout':
      advanced: true

  color_theme:
    type:
    - "null"
    - type: enum
      symbols:
      - "gray"
      - "bw"
      - "linedraw"
      - "light"
      - "dark"
      - "minimal"
      - "classic"
      - "void"
    default: "classic"
    label: "Color theme"
    doc: |
      Color theme for all generated plots.
    'sd:layout':
      advanced: true

  parallel_memory_limit:
    type:
    - "null"
    - type: enum
      symbols:
      - "32"
    default: "32"
    label: "Maximum shared memory in GB"
    doc: |
      Maximum memory in GB allowed to
      be shared between the workers
      when using multiple CPUs.
    'sd:layout':
      advanced: true

  vector_memory_limit:
    type:
    - "null"
    - type: enum
      symbols:
      - "64"
    default: "64"
    label: "Maximum vector memory in GB"
    doc: |
      Maximum vector memory in GB
      allowed to be used by R.
    'sd:layout':
      advanced: true

  threads:
    type:
    - "null"
    - type: enum
      symbols:
      - "1"
    default: "1"
    label: "Number of cores/cpus"
    doc: |
      Number of cores/cpus to use
    'sd:layout':
      advanced: true


outputs:

  umap_rd_rnaumap_plot_png:
    type: File?
    outputSource: de_pseudobulk/umap_rd_rnaumap_plot_png
    label: "Cells RNA UMAP split by selected criteria"
    doc: |
      Cells UMAP split by selected criteria,
      optionally subsetted to the specific
      group (rnaumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Cells RNA UMAP split by selected criteria'

  umap_rd_atacumap_plot_png:
    type: File?
    outputSource: de_pseudobulk/umap_rd_atacumap_plot_png
    label: "Cells ATAC UMAP split by selected criteria"
    doc: |
      Cells UMAP split by selected criteria,
      optionally subsetted to the specific
      group (atacumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Cells ATAC UMAP split by selected criteria'

  umap_rd_wnnumap_plot_png:
    type: File?
    outputSource: de_pseudobulk/umap_rd_wnnumap_plot_png
    label: "Cells WNN UMAP split by selected criteria"
    doc: |
      Cells UMAP split by selected criteria,
      optionally subsetted to the specific
      group (wnnumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Cells WNN UMAP split by selected criteria'

  mds_plot_html:
    type: File?
    outputSource: de_pseudobulk/mds_plot_html
    label: "Interactive MDS Plot"
    doc: |
      MDS plot of pseudobulk aggregated
      normalized reads counts. All genes.
      HTML format
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  volcano_plot_html_file:
    type: File
    outputSource: make_volcano_plot/html_file
    label: "Interactive Volcano Plot"
    doc: |
      HTML index file for Volcano Plot
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  volcano_plot_html_data:
    type: Directory
    outputSource: make_volcano_plot/html_data
    label: "Directory html data for Volcano Plot"
    doc: |
      Directory html data for Volcano Plot

  heatmap_html:
    type: File
    outputSource: morpheus_heatmap/heatmap_html
    label: "Interactive Gene Expression Heatmap"
    doc: |
      Morpheus heatmap in HTML format
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  pca_1_2_plot_png:
    type: File?
    outputSource: de_pseudobulk/pca_1_2_plot_png
    label: "Normalized reads counts PCA (1, 2). All genes."
    doc: |
      Normalized reads counts PCA (1, 2). All genes.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Normalized reads counts PCA (1, 2). All genes'
  
  pca_2_3_plot_png:
    type: File?
    outputSource: de_pseudobulk/pca_2_3_plot_png
    label: "Normalized reads counts PCA (2, 3). All genes."
    doc: |
      Normalized reads counts PCA (2, 3). All genes.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Normalized reads counts PCA (2, 3). All genes'

  dxpr_vlcn_plot_png:
    type: File?
    outputSource: de_pseudobulk/dxpr_vlcn_plot_png
    label: "Volcano plot of differentially expressed genes"
    doc: |
      Volcano plot of differentially expressed genes.
      Highlighed genes are either provided by user or
      top 10 genes with the highest log2FoldChange
      values. The direction of comparison is defined
      as --second vs --first. Cells are optionally
      subsetted to the specific group and optionally
      coerced to the pseudobulk form.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Volcano plot of differentially expressed genes'

  xpr_dnst_plot_png:
    type: File?
    outputSource: de_pseudobulk/xpr_dnst_plot_png
    label: "Log normalized gene expression density plots"
    doc: |
      Log normalized gene expression density plots for
      either user provided or top 10 differentially
      expressed genes with the highest log2FoldChange
      values. The direction of comparison is defined
      as --second vs --first.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Log normalized gene expression density plots'

  xpr_per_cell_rd_rnaumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: de_pseudobulk/xpr_per_cell_rd_rnaumap_plot_png
    label: "Log normalized gene expression on cells RNA UMAP"
    doc: |
      Log normalized gene expression on cells UMAP
      split by selected criteria, optionally subsetted
      to the specific group (rnaumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression RNA'
        Caption: 'Log normalized gene expression on cells RNA UMAP'

  xpr_per_cell_rd_atacumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: de_pseudobulk/xpr_per_cell_rd_atacumap_plot_png
    label: "Log normalized gene expression on cells ATAC UMAP"
    doc: |
      Log normalized gene expression on cells UMAP
      split by selected criteria, optionally subsetted
      to the specific group (atacumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression ATAC'
        Caption: 'Log normalized gene expression on cells ATAC UMAP'

  xpr_per_cell_rd_wnnumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: de_pseudobulk/xpr_per_cell_rd_wnnumap_plot_png
    label: "Log normalized gene expression on cells WNN UMAP"
    doc: |
      Log normalized gene expression on cells UMAP
      split by selected criteria, optionally subsetted
      to the specific group (wnnumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression WNN'
        Caption: 'Log normalized gene expression on cells WNN UMAP'

  xpr_htmp_plot_png:
    type: File?
    outputSource: de_pseudobulk/xpr_htmp_plot_png
    label: "Filtered by adjusted P-value normalized gene expression heatmap"
    doc: |
      Filtered by adjusted P-value normalized gene
      expression heatmap per cell optionally subsetted
      to the specific group.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Filtered by adjusted P-value normalized gene expression heatmap'

  diff_expr_genes:
    type: File
    outputSource: de_pseudobulk/diff_expr_genes
    label: "Differentially expressed genes. Not filtered"
    doc: |
      Differentially expressed genes. Not filtered
      by adjusted P-value.
      TSV format
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Diff expressed genes'
        Title: 'Differentially expressed genes. Not filtered'

  read_counts_file:
    type: File?
    outputSource: de_pseudobulk/bulk_read_counts_gct
    label: "GSEA compatible not filtered normalized reads counts"
    doc: |
      GSEA compatible not filtered normalized reads
      counts aggregated to pseudobulk form.
      GCT format

  phenotypes_file:
    type: File?
    outputSource: de_pseudobulk/bulk_phenotypes_cls
    label: "GSEA compatible phenotypes file"
    doc: |
      GSEA compatible phenotypes file defined based
      on --splitby, --first, and --second parameters.
      CLS format

  cell_read_counts_gct:
    type: File
    outputSource: de_pseudobulk/cell_read_counts_gct
    label: "Filtered normalized reads counts per cell"
    doc: |
      Filtered normalized reads counts per cell.
      GCT format

  de_pseudobulk_stdout_log:
    type: File
    outputSource: de_pseudobulk/stdout_log
    label: "stdout log generated by de_pseudobulk step"
    doc: |
      stdout log generated by de_pseudobulk step

  de_pseudobulk_stderr_log:
    type: File
    outputSource: de_pseudobulk/stderr_log
    label: "stderr log generated by de_pseudobulk step"
    doc: |
      stderr log generated by de_pseudobulk step

  morpheus_heatmap_stdout_log:
    type: File
    outputSource: morpheus_heatmap/stdout_log
    label: "stdout log generated by morpheus_heatmap step"
    doc: "stdout log generated by morpheus_heatmap step"

  morpheus_heatmap_stderr_log:
    type: File
    outputSource: morpheus_heatmap/stderr_log
    label: "stderr log generated by morpheus_heatmap step"
    doc: "stderr log generated by morpheus_heatmap step"


steps:

  de_pseudobulk:
    run: ../tools/sc-rna-de-pseudobulk.cwl
    in:
      query_data_rds: query_data_rds
      datasets_metadata: datasets_metadata
      barcodes_data: barcodes_data
      groupby:
        source: groupby
        valueFrom: $(self==""?null:self)            # safety measure
      subset:
        source: subset
        valueFrom: $(split_features(self))
      splitby: splitby
      first_cond: first_cond
      second_cond: second_cond
      analysis_method:
        source: analysis_method
        valueFrom: $(self.split(" ")[0])
      batchby:
        source: batchby
        valueFrom: $(self==""?null:self)            # safety measure
      maximum_padj: maximum_padj
      genes_of_interest:
        source: genes_of_interest
        valueFrom: $(split_features(self))
      exclude_pattern:
        source: exclude_pattern
        valueFrom: $(self==""?null:self)            # safety measure
      cluster_method:
        source: cluster_method
        valueFrom: $(self=="none"?null:self)
      row_distance: row_distance
      column_distance: column_distance
      center_row: center_row
      color_theme: color_theme
      verbose:
        default: true
      parallel_memory_limit:
        source: parallel_memory_limit
        valueFrom: $(parseInt(self))
      vector_memory_limit:
        source: vector_memory_limit
        valueFrom: $(parseInt(self))
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
      - umap_rd_rnaumap_plot_png
      - umap_rd_atacumap_plot_png
      - umap_rd_wnnumap_plot_png
      - mds_plot_html
      - pca_1_2_plot_png
      - pca_2_3_plot_png
      - dxpr_vlcn_plot_png
      - xpr_dnst_plot_png
      - xpr_per_cell_rd_rnaumap_plot_png
      - xpr_per_cell_rd_atacumap_plot_png
      - xpr_per_cell_rd_wnnumap_plot_png
      - xpr_htmp_plot_png
      - diff_expr_genes
      - bulk_read_counts_gct
      - bulk_phenotypes_cls
      - cell_read_counts_gct
      - stdout_log
      - stderr_log

  morpheus_heatmap:
    run: ../tools/morpheus-heatmap.cwl
    in:
     read_counts_gct: de_pseudobulk/cell_read_counts_gct
    out:
    - heatmap_html
    - stdout_log
    - stderr_log

  make_volcano_plot:
    run: ../tools/volcano-plot.cwl
    in:
      diff_expr_file: de_pseudobulk/diff_expr_genes
      x_axis_column:
        default: "log2FoldChange"
      y_axis_column:
        default: "padj"
      label_column:
        default: "gene"
    out:
    - html_data
    - html_file


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell Pseudobulk Differential Expression Analysis Between Datasets"
s:name: "Single-cell Pseudobulk Differential Expression Analysis Between Datasets"
s:alternateName: "Identifies differentially expressed genes between groups of cells coerced to pseudobulk datasets"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-rna-de-pseudobulk.cwl
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
  Single-cell Pseudobulk Differential Expression Analysis Between Datasets

  Identifies differentially expressed genes between groups of cells
  coerced to pseudobulk datasets.