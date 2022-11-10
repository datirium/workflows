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


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Experiment run through any of the Single-cell Cluster or Manual Cell Type Assignment Analysis"
    doc: |
      Path to the RDS file to load Seurat object from. This file should include genes
      expression information stored in the RNA assay. Additionally, 'rnaumap', and/or
      'atacumap', and/or 'wnnumap' dimensionality reductions should be present.
    'sd:upstreamSource': "sc_tools_sample/seurat_data_rds"
    'sd:localLabel': true

  splitby:
    type: string
    label: "Column from the Seurat object metadata to split datasets into two groups"
    doc: |
      Column from the Seurat object metadata to split datasets into two groups
      to run --second vs --first pseudobulk DE analysis, i.e., calculate log2FC.
      May be one of the columns from the extra metadata added with --metadata
      parameter. Provided value should group the datasets, not cells, therefore
      do not use a column with clustering results.

  first_cond:
    type: string
    label: "Value from the Seurat object metadata column to define the first group of datasets"
    doc: |
      Value from the Seurat object metadata column set with --splitby to define the
      first group of datasets for pseudobulk DE analysis.

  second_cond:
    type: string
    label: "Value from the Seurat object metadata column to define the second group of datasets"
    doc: |
      Value from the Seurat object metadata column set with --splitby to define the
      second group of datasets for pseudobulk DE analysis.

  batchby:
    type: string?
    default: null
    label: "Column from the Seurat object metadata to group datasets into batches"
    doc: |
      Column from the Seurat object metadata to group datasets into batches. It will be used
      as a factor variable to model batch effect when running pseudobulk DE analysis (makes
      design formula look like ~splitby+batchby). May be one of the columns from the extra
      metadata added with --metadata parameter. Provided value should batch the datasets, not
      cells, therefore do not use a column with clustering results. Default: do not model
      batch effect.

  groupby:
    type: string?
    default: null
    label: "Column from the Seurat object metadata to group cells for optional subsetting"
    doc: |
      Column from the Seurat object metadata to group cells for optional subsetting
      when combined with --subset parameter. May be one of the columns from the extra
      metadata added with --metadata parameter. Ignored if --subset is not set. Provided
      value defines the groups of cells, therefore any metadata column, including the
      clustering results, may be used. Default: do not subset, run pseudobulk DE analysis
      for all cells jointly

  subset:
    type: string?
    default: null
    label: "Value(s) to subset cells before running analysis"
    doc: |
      Value(s) from the column set with --groupby parameter to subset cells
      before running pseudobulk DE analysis. If multiple values are provided
      run analysis jointly for selected groups of cells. Ignored if --groupby
      is not set. Default: do not subset, run pseudobulk DE analysis for all
      cells jointly

  datasets_metadata:
    type: File?
    label: "Path to the TSV/CSV file to optionally extend Seurat object metadata"
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat object metadata with
      categorical values using samples identities. First column - 'library_id'
      should correspond to all unique values from the 'new.ident' column of the
      loaded Seurat object. If any of the provided in this file columns are already
      present in the Seurat object metadata, they will be overwritten. Default: no
      extra metadata is added

  lrt:
    type: boolean?
    default: false
    label: "Use LRT instead of the pair-wise Wald test"
    doc: |
      Use LRT instead of the pair-wise Wald test. If --batchby is not provided
      use ~1 as a reduced formula, otherwise ~batchby. Default: use Wald test
    'sd:layout':
      advanced: true

  maximum_padj:
    type: float?
    default: 0.05
    label: "Maximum significance level used in the exploratory visualization part of the analysis"
    doc: |
      In the exploratory visualization part of the analysis output only features
      with adjusted P-value not bigger than this value. Default: 0.05
    'sd:layout':
      advanced: true

  genes_of_interest:
    type: string?
    default: null
    label: "Genes of interest to label on the generated plots"
    doc: |
      Genes of interest to label on the generated plots. Default: top 10 genes
      with the highest and the lowest log2FC expression values.
    'sd:layout':
      advanced: true

  exclude_pattern:
    type: string?
    default: null
    label: "Regex pattern to identify and exclude non-coding RNA genes from the analysis"
    doc: |
      Regex pattern to identify and exclude non-coding RNA genes from the pseudobulk
      DE analysis (not case-sensitive). If any of such genes were provided in the --genes
      parameter, they will be excluded from there as well.
    'sd:layout':
      advanced: true

  normalization_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "vst"
      - "rlog"
    default: "rlog"
    label: "Read counts normalization for the exploratory visualization part of the analysis"
    doc: |
      Read counts normalization for the exploratory visualization part of the analysis.
      Use 'vst' for medium-to-large datasets (n > 30) and 'rlog' for small datasets
      (n < 30), when there is a wide range of sequencing depth across samples.
      Default: rlog
    'sd:layout':
      advanced: true

  remove:
    type: boolean?
    default: false
    label: "Remove batch effect when generating normalized read counts"
    doc: |
      Remove batch effect when generating normalized read counts for the exploratory
      visualization part of the analysis. Ignored if --batchby is not provided.
      Default: do not remove batch effect from normalized read counts.
    'sd:layout':
      advanced: true

  center_row:
    type: boolean?
    default: false
    label: "Apply mean centering for feature expression prior to running clustering by row"
    doc: |
      Apply mean centering for gene expression prior to running
      clustering by row. Ignored if --cluster is set to column or
      not provided. Default: do not centered
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
    default: "none"
    label: "Hopach clustering method to be run on normalized read counts"
    doc: |
      Hopach clustering method to be run on normalized read counts for the
      exploratory visualization part of the analysis. Default: do not run
      clustering
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
    label: "Distance metric for HOPACH row clustering"
    doc: |
      Distance metric for HOPACH row clustering. Ignored if --cluster is set
      to column or not provided. Default: cosangle
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
    label: "Distance metric for HOPACH column clustering"
    doc: |
      Distance metric for HOPACH column clustering. Ignored if --cluster is set
      to row or not provided. Default: euclid
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
    label: "Color theme for all generated plots"
    doc: |
      Color theme for all generated plots. One of gray, bw, linedraw, light,
      dark, minimal, classic, void.
      Default: classic
    'sd:layout':
      advanced: true

  parallel_memory_limit:
    type:
    - "null"
    - type: enum
      symbols:
      - "32"
    default: "32"
    label: "Maximum memory in GB allowed to be shared between the workers when using multiple CPUs"
    doc: |
      Maximum memory in GB allowed to be shared between the workers
      when using multiple --cpus.
      Forced to 32 GB
    'sd:layout':
      advanced: true

  vector_memory_limit:
    type:
    - "null"
    - type: enum
      symbols:
      - "64"
    default: "64"
    label: "Maximum vector memory in GB allowed to be used by R"
    doc: |
      Maximum vector memory in GB allowed to be used by R.
      Forced to 64 GB
    'sd:layout':
      advanced: true

  threads:
    type:
    - "null"
    - type: enum
      symbols:
      - "1"
    default: "1"
    label: "Number of cores/cpus to use"
    doc: |
      Number of cores/cpus to use
      Forced to 1
    'sd:layout':
      advanced: true


outputs:

  umap_rd_rnaumap_plot_png:
    type: File?
    outputSource: de_pseudobulk/umap_rd_rnaumap_plot_png
    label: "Cells RNA UMAP split by selected biological condition"
    doc: |
      Cells UMAP split by selected biological condition, optionally
      subsetted to the specific cluster or cell type (rnaumap dim.
      reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Cells RNA UMAP split by selected biological condition'

  umap_rd_atacumap_plot_png:
    type: File?
    outputSource: de_pseudobulk/umap_rd_atacumap_plot_png
    label: "Cells ATAC UMAP split by selected biological condition"
    doc: |
      Cells UMAP split by selected biological condition, optionally
      subsetted to the specific cluster or cell type (atacumap dim.
      reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Cells ATAC UMAP split by selected biological condition'

  umap_rd_wnnumap_plot_png:
    type: File?
    outputSource: de_pseudobulk/umap_rd_wnnumap_plot_png
    label: "Cells WNN UMAP split by selected biological condition"
    doc: |
      Cells UMAP split by selected biological condition, optionally
      subsetted to the specific cluster or cell type (wnnumap dim.
      reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Cells WNN UMAP split by selected biological condition'

  mds_plot_html:
    type: File?
    outputSource: de_pseudobulk/mds_plot_html
    label: "MDS plot of normalized counts"
    doc: |
      MDS plot of normalized counts. Optionally batch corrected
      if --remove was set to True.
      HTML format
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  heatmap_html:
    type: File
    outputSource: morpheus_heatmap/heatmap_html
    label: "Heatmap of normalized counts"
    doc: |
      Morpheus heatmap in HTML format
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  pca_1_2_plot_png:
    type: File?
    outputSource: de_pseudobulk/pca_1_2_plot_png
    label: "Normalized counts PCA (PC1 and PC2)"
    doc: |
      Normalized counts PCA (PC1 and PC2) subsetted to all DE genes regardless
      of Padj, optionally batch corrected by the selected criteria.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Normalized counts PCA (PC1 and PC2)'
  
  pca_2_3_plot_png:
    type: File?
    outputSource: de_pseudobulk/pca_2_3_plot_png
    label: "Normalized counts PCA (PC2 and PC3)"
    doc: |
      Normalized counts PCA (PC2 and PC3) subsetted to all DE genes regardless
      of Padj, optionally batch corrected by the selected criteria.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Normalized counts PCA (PC2 and PC3)'

  dxpr_vlcn_plot_png:
    type: File?
    outputSource: de_pseudobulk/dxpr_vlcn_plot_png
    label: "Volcano plot of differentially expressed genes"
    doc: |
      Volcano plot of differentially expressed genes. Highlighed genes are either
      provided by user or top 10 genes with the highest log2FC values. The direction
      of comparison is defined by --second vs --first groups of cells optionally
      subsetted to the specific cluster or cell type and coerced to the pseudobulk
      RNA-Seq samples.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Volcano plot of differentially expressed genes'

  xpr_dnst_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: de_pseudobulk/xpr_dnst_plot_png
    label: "Log normalized gene expression density per dataset"
    doc: |
      Log normalized gene expression density per dataset optionally subsetted to the
      specific cluster or cell type.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Log normalized gene expression density per dataset'

  xpr_per_cell_rd_rnaumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: de_pseudobulk/xpr_per_cell_rd_rnaumap_plot_png
    label: "Log normalized gene expression on cells RNA UMAP per dataset"
    doc: |
      Log normalized gene expression on cells UMAP per dataset optionally subsetted
      to the specific cluster or cell type (rnaumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Log normalized gene expression on cells RNA UMAP per dataset'

  xpr_per_cell_rd_atacumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: de_pseudobulk/xpr_per_cell_rd_atacumap_plot_png
    label: "Log normalized gene expression on cells ATAC UMAP per dataset"
    doc: |
      Log normalized gene expression on cells UMAP per dataset optionally subsetted
      to the specific cluster or cell type (atacumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Log normalized gene expression on cells ATAC UMAP per dataset'

  xpr_per_cell_rd_wnnumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: de_pseudobulk/xpr_per_cell_rd_wnnumap_plot_png
    label: "Log normalized gene expression on cells WNN UMAP per dataset"
    doc: |
      Log normalized gene expression on cells UMAP per dataset optionally subsetted
      to the specific cluster or cell type (wnnumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Log normalized gene expression on cells WNN UMAP per dataset'

  xpr_htmp_plot_png:
    type: File?
    outputSource: de_pseudobulk/xpr_htmp_plot_png
    label: "Log normalized gene expression heatmap per dataset"
    doc: |
      Normalized gene expression heatmap optionally subsetted
      to the specific cluster or cell type.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Normalized gene expression heatmap'

  diff_expr_genes:
    type: File
    outputSource: de_pseudobulk/diff_expr_genes
    label: "Differentially expressed genes"
    doc: |
      Differentially expressed genes.
      TSV format
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Diff expressed genes'
        Title: 'Differentially expressed genes'

  read_counts_file:
    type: File
    outputSource: de_pseudobulk/read_counts_gct
    label: "GSEA compatible normalized counts"
    doc: |
      GSEA compatible normalized counts, optionally, batch corrected.
      GCT format

  phenotypes_file:
    type: File
    outputSource: de_pseudobulk/phenotypes_cls
    label: "GSEA compatible phenotypes file"
    doc: |
      GSEA compatible phenotypes file defined based on --splitby, --first,
      and --second parameters.
      CLS format

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
      splitby: splitby
      first_cond: first_cond
      second_cond: second_cond
      batchby: batchby
      groupby: groupby
      subset:
        source: subset
        valueFrom: $(split_features(self))
      lrt: lrt
      maximum_padj: maximum_padj
      genes_of_interest:
        source: genes_of_interest
        valueFrom: $(split_features(self))
      exclude_pattern: exclude_pattern
      normalization_method: normalization_method
      remove: remove
      cluster_method:
        source: cluster_method
        valueFrom: $(self=="none"?null:self)
      row_distance: row_distance
      column_distance: column_distance
      center_row: center_row
      verbose:
        default: true
      color_theme: color_theme
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
      - read_counts_gct
      - phenotypes_cls
      - stdout_log
      - stderr_log

  morpheus_heatmap:
    run: ../tools/morpheus-heatmap.cwl
    in:
     read_counts_gct: de_pseudobulk/read_counts_gct
    out:
    - heatmap_html
    - stdout_log
    - stderr_log


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