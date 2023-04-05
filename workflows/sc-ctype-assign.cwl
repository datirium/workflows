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
    - var get_query_column = function(prefix, reduction, resolution) {
          if (reduction=="RNA") {
            return prefix + "rna_res." + resolution;
          } else if (reduction=="ATAC") {
            return prefix + "atac_res." + resolution;
          } else if (reduction=="WNN") {
            return prefix + "wsnn_res." + resolution;
          }
      };


'sd:upstream':
  sc_tools_sample:
  - "sc-ctype-assign.cwl"
  - "sc-rna-cluster.cwl"
  - "sc-atac-cluster.cwl"
  - "sc-wnn-cluster.cwl"
  sc_arc_sample:
  - "cellranger-arc-count.cwl"
  - "cellranger-arc-aggr.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Experiment run through any of the Single-cell Cluster Analysis"
    doc: |
      Path to the RDS file to load Seurat object from. This file should include
      genes expression and/or chromatin accessibility information stored in the RNA
      and ATAC assays correspondingly. Additionally, 'rnaumap', and/or 'atacumap',
      and/or 'wnnumap' dimensionality reductions should be present.
    'sd:upstreamSource': "sc_tools_sample/seurat_data_rds"
    'sd:localLabel': true

  query_reduction:
    type:
    - "null"
    - type: enum
      symbols:
      - "RNA"
      - "ATAC"
      - "WNN"
    default: "RNA"
    label: "Select clusters based on"
    doc: |
      If set to 'RNA', then 'get_query_column' will have suffix 'rna_res'.
      If set to 'ATAC', then 'get_query_column' will have suffix 'atac_res'.
      If set to 'WNN', then 'get_query_column' will have suffix 'wsnn_res'.
  
  query_resolution:
    type: float
    label: "Clustering resolution to assign cell types to"
    doc: |
      Clustering resolution defines 'query_source_column' and 'query_target_column'
      inputs for 'assign_cell_types' step

  atac_fragments_file:
    type: File?
    secondaryFiles:
    - .tbi
    label: "Cell Ranger ARC Count/Aggregate Experiment for ATAC or WNN clusters"
    doc: |
      Count and barcode information for every ATAC fragment used in the loaded Seurat
      object. File should be saved in TSV format with tbi-index file. Ignored if the
      loaded Seurat object doesn't include ATAC assay.
    'sd:upstreamSource': "sc_arc_sample/atac_fragments_file"
    'sd:localLabel': true

  genes_of_interest:
    type: string?
    default: null
    label: "Genes of interest to build gene expression and/or Tn5 insertion frequency plots for the nearest peaks"
    doc: |
      Genes of interest to build gene expression and/or Tn5 insertion frequency plots
      for the nearest peaks. To build gene expression plots the loaded Seurat object
      should include RNA assay. To build Tn5 insertion frequency plots for the nearest
      peaks the loaded Seurat object should include ATAC assay as well as the --fragments
      file should be provided.
      Default: None

  cell_type_data:
    type: File
    label: "TSV/CSV cell types metadata file with 'cluster' and 'type' columns"
    doc: |
      Path to the TSV/CSV file for manual cell type assignment for each of the clusters.
      First column - 'cluster', second column may have arbitrary name.

  identify_diff_genes:
    type: boolean?
    default: false
    label: "Identify differentially expressed genes for assigned cell types"
    doc: |
      Identify differentially expressed genes (putative gene markers) for
      assigned cell types. Ignored if loaded Seurat object doesn't include
      genes expression information stored in the RNA assay.
      Default: false
    'sd:layout':
      advanced: true

  identify_diff_peaks:
    type: boolean?
    default: false
    label: "Identify differentially accessible peaks for assigned cell types"
    doc: |
      Identify differentially accessible peaks for assigned cell types. Ignored
      if loaded Seurat object doesn't include chromatin accessibility information
      stored in the ATAC assay.
      Default: false
    'sd:layout':
      advanced: true

  rna_minimum_logfc:
    type: float?
    default: 0.25
    label: "Include only those genes that on average have log fold change difference in expression between every tested pair of cell types not lower than this value"
    doc: |
      For putative gene markers identification include only those genes that
      on average have log fold change difference in expression between every
      tested pair of cell types not lower than this value. Ignored if '--diffgenes'
      is not set or RNA assay is not present.
      Default: 0.25
    'sd:layout':
      advanced: true

  rna_minimum_pct:
    type: float?
    default: 0.1
    label: "Include only those genes that are detected in not lower than this fraction of cells in either of the two tested cell types"
    doc: |
      For putative gene markers identification include only those genes that
      are detected in not lower than this fraction of cells in either of the
      two tested cell types. Ignored if '--diffgenes' is not set or RNA assay
      is not present.
      Default: 0.1
    'sd:layout':
      advanced: true

  atac_minimum_logfc:
    type: float?
    default: 0.25
    label: "Include only those peaks that on average have log fold change difference in the chromatin accessibility between every tested pair of cell types not lower than this value"
    doc: |
      For differentially accessible peaks identification include only those peaks that
      on average have log fold change difference in the chromatin accessibility between
      every tested pair of cell types not lower than this value. Ignored if '--diffpeaks'
      is not set or ATAC assay is not present.
      Default: 0.25
    'sd:layout':
      advanced: true

  atac_minimum_pct:
    type: float?
    default: 0.05
    label: "Include only those peaks that are detected in not lower than this fraction of cells in either of the two tested cell types"
    doc: |
      For differentially accessible peaks identification include only those peaks that
      are detected in not lower than this fraction of cells in either of the two tested
      cell types. Ignored if '--diffpeaks' is not set or ATAC assay is not present.
      Default: 0.05
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
    outputSource: ctype_assign/umap_rd_rnaumap_plot_png
    label: "Clustered cells RNA UMAP with assigned cell types"
    doc: |
      Cells UMAP with assigned cell types (rnaumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Clustered cells RNA UMAP with assigned cell types'

  umap_rd_atacumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_rd_atacumap_plot_png
    label: "Clustered cells ATAC UMAP with assigned cell types"
    doc: |
      Cells UMAP with assigned cell types (atacumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Clustered cells ATAC UMAP with assigned cell types'

  umap_rd_wnnumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_rd_wnnumap_plot_png
    label: "Clustered cells WNN UMAP with assigned cell types"
    doc: |
      Cells UMAP with assigned cell types (wnnumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Clustered cells WNN UMAP with assigned cell types'

  umap_spl_idnt_rd_rnaumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_idnt_rd_rnaumap_plot_png
    label: "Split by dataset clustered cells RNA UMAP with assigned cell types"
    doc: |
      Split by dataset cells UMAP with assigned cell types (rnaumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Split by dataset clustered cells RNA UMAP with assigned cell types'

  umap_spl_idnt_rd_atacumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_idnt_rd_atacumap_plot_png
    label: "Split by dataset clustered cells ATAC UMAP with assigned cell types"
    doc: |
      Split by dataset cells UMAP with assigned cell types (atacumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Split by dataset clustered cells ATAC UMAP with assigned cell types'

  umap_spl_idnt_rd_wnnumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_idnt_rd_wnnumap_plot_png
    label: "Split by dataset clustered cells WNN UMAP with assigned cell types"
    doc: |
      Split by dataset cells UMAP with assigned cell types (wnnumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Split by dataset clustered cells WNN UMAP with assigned cell types'

  umap_spl_cnd_rd_rnaumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_cnd_rd_rnaumap_plot_png
    label: "Split by grouping condition clustered cells RNA UMAP with assigned cell types"
    doc: |
      Split by grouping condition cells UMAP with assigned cell types (rnaumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Split by grouping condition clustered cells RNA UMAP with assigned cell types'

  umap_spl_cnd_rd_atacumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_cnd_rd_atacumap_plot_png
    label: "Split by grouping condition clustered cells ATAC UMAP with assigned cell types"
    doc: |
      Split by grouping condition cells UMAP with assigned cell types (atacumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Split by grouping condition clustered cells ATAC UMAP with assigned cell types'

  umap_spl_cnd_rd_wnnumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_cnd_rd_wnnumap_plot_png
    label: "Split by grouping condition clustered cells WNN UMAP with assigned cell types"
    doc: |
      Split by grouping condition cells UMAP with assigned cell types (wnnumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Split by grouping condition clustered cells WNN UMAP with assigned cell types'

  umap_spl_ph_rd_rnaumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_ph_rd_rnaumap_plot_png
    label: "Split by cell cycle phase cells RNA UMAP with assigned cell types"
    doc: |
      Split by cell cycle phase cells UMAP with assigned cell types (rnaumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Split by cell cycle phase cells RNA UMAP with assigned cell types'

  umap_spl_ph_rd_atacumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_ph_rd_atacumap_plot_png
    label: "Split by cell cycle phase cells ATAC UMAP with assigned cell types"
    doc: |
      Split by cell cycle phase cells UMAP with assigned cell types (atacumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Split by cell cycle phase cells ATAC UMAP with assigned cell types'

  umap_spl_ph_rd_wnnumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_ph_rd_wnnumap_plot_png
    label: "Split by cell cycle phase cells WNN UMAP with assigned cell types"
    doc: |
      Split by cell cycle phase cells UMAP with assigned cell types (wnnumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Split by cell cycle phase cells WNN UMAP with assigned cell types'

  cmp_gr_ctyp_spl_idnt_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_ctyp_spl_idnt_plot_png
    label: "Grouped by cell type split by dataset cells composition plot. Downsampled."
    doc: |
      Grouped by cell type split by dataset cells composition plot. Downsampled.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Grouped by cell type split by dataset cells composition plot. Downsampled.'

  cmp_gr_idnt_spl_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_idnt_spl_ctyp_plot_png
    label: "Grouped by dataset split by cell type cells composition plot. Downsampled."
    doc: |
      Grouped by dataset split by cell type cells composition plot. Downsampled.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Grouped by dataset split by cell type cells composition plot. Downsampled.'

  cmp_gr_ph_spl_idnt_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_ph_spl_idnt_plot_png
    label: "Grouped by cell cycle phase split by dataset cells composition plot. Downsampled."
    doc: |
      Grouped by cell cycle phase split by dataset cells composition plot. Downsampled.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Grouped by cell cycle phase split by dataset cells composition plot. Downsampled.'

  cmp_gr_ctyp_spl_cnd_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_ctyp_spl_cnd_plot_png
    label: "Grouped by cell type split by condition cells composition plot. Downsampled."
    doc: |
      Grouped by cell type split by condition cells composition plot. Downsampled.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Grouped by cell type split by condition cells composition plot. Downsampled.'

  cmp_gr_cnd_spl_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_cnd_spl_ctyp_plot_png
    label: "Grouped by condition split by cell type cells composition plot. Downsampled."
    doc: |
      Grouped by condition split by cell type cells composition plot. Downsampled.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Grouped by condition split by cell type cells composition plot. Downsampled.'

  cmp_gr_ph_spl_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_ph_spl_ctyp_plot_png
    label: "Grouped by cell cycle phase split by cell type cells composition plot. Downsampled."
    doc: |
      Grouped by cell cycle phase split by cell type cells composition plot. Downsampled.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Grouped by cell cycle phase split by cell type cells composition plot. Downsampled.'

  xpr_avg_plot_png:
    type: File?
    outputSource: ctype_assign/xpr_avg_plot_png
    label: "Log normalized scaled average gene expression per cell type"
    doc: |
      Log normalized scaled average gene expression per cell type.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Log normalized scaled average gene expression per cell type'

  xpr_dnst_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: ctype_assign/xpr_dnst_plot_png
    label: "Log normalized gene expression density per cell type"
    doc: |
      Log normalized gene expression density per cell type.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Log normalized gene expression density per cell type'

  xpr_per_cell_rd_rnaumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: ctype_assign/xpr_per_cell_rd_rnaumap_plot_png
    label: "Log normalized gene expression on cells RNA UMAP with assigned cell types"
    doc: |
      Log normalized gene expression on cells UMAP with assigned cell types (rnaumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Log normalized gene expression on cells RNA UMAP with assigned cell types'

  xpr_per_cell_rd_atacumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: ctype_assign/xpr_per_cell_rd_atacumap_plot_png
    label: "Log normalized gene expression on cells ATAC UMAP with assigned cell types"
    doc: |
      Log normalized gene expression on cells UMAP with assigned cell types (atacumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Log normalized gene expression on cells ATAC UMAP with assigned cell types'

  xpr_per_cell_rd_wnnumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: ctype_assign/xpr_per_cell_rd_wnnumap_plot_png
    label: "Log normalized gene expression on cells WNN UMAP with assigned cell types"
    doc: |
      Log normalized gene expression on cells UMAP with assigned cell types (wnnumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Log normalized gene expression on cells WNN UMAP with assigned cell types'

  xpr_per_cell_sgnl_rd_rnaumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: ctype_assign/xpr_per_cell_sgnl_rd_rnaumap_plot_png
    label: "Log normalized gene expression density on cells RNA UMAP with assigned cell types"
    doc: |
      Log normalized gene expression density on cells UMAP with assigned cell types (rnaumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Log normalized gene expression density on cells RNA UMAP with assigned cell types'

  xpr_per_cell_sgnl_rd_atacumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: ctype_assign/xpr_per_cell_sgnl_rd_atacumap_plot_png
    label: "Log normalized gene expression density on cells ATAC UMAP with assigned cell types"
    doc: |
      Log normalized gene expression density on cells UMAP with assigned cell types (atacumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Log normalized gene expression density on cells ATAC UMAP with assigned cell types'

  xpr_per_cell_sgnl_rd_wnnumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: ctype_assign/xpr_per_cell_sgnl_rd_wnnumap_plot_png
    label: "Log normalized gene expression density on cells WNN UMAP with assigned cell types"
    doc: |
      Log normalized gene expression density on cells UMAP with assigned cell types (wnnumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Log normalized gene expression density on cells WNN UMAP with assigned cell types'

  cvrg_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: ctype_assign/cvrg_plot_png
    label: "Tn5 insertion frequency plot around gene"
    doc: |
      Tn5 insertion frequency plot around gene.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Genome coverage'
        Caption: 'Tn5 insertion frequency plot around gene'

  xpr_htmp_plot_png:
    type: File?
    outputSource: ctype_assign/xpr_htmp_plot_png
    label: "Normalized gene expression heatmap grouped by cell type"
    doc: |
      Normalized gene expression heatmap grouped by cell type.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression'
        Caption: 'Normalized gene expression heatmap grouped by cell type'

  gene_markers_tsv:
    type: File?
    outputSource: ctype_assign/gene_markers_tsv
    label: "Differentially expressed genes between each pair of cell types"
    doc: |
      Differentially expressed genes between each pair of cell types.
      TSV format
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Gene markers'
        Title: 'Differentially expressed genes between each pair of cell types'

  peak_markers_tsv:
    type: File?
    outputSource: ctype_assign/peak_markers_tsv
    label: "Differentially accessible peaks between each pair of cell types"
    doc: |
      Differentially accessible peaks between each pair of cell types.
      TSV format
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Diff. peaks'
        Title: 'Differentially accessible peaks between each pair of cell types'

  ucsc_cb_config_data:
    type: File
    outputSource: compress_cellbrowser_config_data/compressed_folder
    label: "Compressed directory with UCSC Cellbrowser configuration data"
    doc: |
      Compressed directory with UCSC Cellbrowser configuration data.

  ucsc_cb_html_data:
    type: Directory
    outputSource: ctype_assign/ucsc_cb_html_data
    label: "Directory with UCSC Cellbrowser html data"
    doc: |
      Directory with UCSC Cellbrowser html data.

  ucsc_cb_html_file:
    type: File
    outputSource: ctype_assign/ucsc_cb_html_file
    label: "Open in UCSC Cell Browser"
    doc: |
      HTML index file from the directory with UCSC Cellbrowser html data.
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  seurat_data_rds:
    type: File
    outputSource: ctype_assign/seurat_data_rds
    label: "Processed Seurat data in RDS format"
    doc: |
      Processed Seurat data in RDS format

  ctype_assign_stdout_log:
    type: File
    outputSource: ctype_assign/stdout_log
    label: "stdout log generated by ctype_assign step"
    doc: |
      stdout log generated by ctype_assign step

  ctype_assign_stderr_log:
    type: File
    outputSource: ctype_assign/stderr_log
    label: "stderr log generated by ctype_assign step"
    doc: |
      stderr log generated by ctype_assign step


steps:

  ctype_assign:
    run: ../tools/sc-ctype-assign.cwl
    in:
      query_data_rds: query_data_rds
      cell_type_data: cell_type_data
      query_source_column:
        source: [query_reduction, query_resolution]
        valueFrom: $(get_query_column("", self[0], self[1]))
      query_target_column:
        source: [query_reduction, query_resolution]
        valueFrom: $(get_query_column("custom_", self[0], self[1]))
      atac_fragments_file: atac_fragments_file
      genes_of_interest:
        source: genes_of_interest
        valueFrom: $(split_features(self))
      identify_diff_genes: identify_diff_genes
      identify_diff_peaks: identify_diff_peaks
      rna_minimum_logfc: rna_minimum_logfc
      rna_minimum_pct: rna_minimum_pct
      atac_minimum_logfc: atac_minimum_logfc
      atac_minimum_pct: atac_minimum_pct
      only_positive_diff_genes:
        default: true
      rna_test_to_use: 
        default: wilcox
      atac_test_to_use:
        default: LR
      verbose:
        default: true
      export_ucsc_cb:
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
    - umap_spl_idnt_rd_rnaumap_plot_png
    - umap_spl_idnt_rd_atacumap_plot_png
    - umap_spl_idnt_rd_wnnumap_plot_png
    - umap_spl_cnd_rd_rnaumap_plot_png
    - umap_spl_cnd_rd_atacumap_plot_png
    - umap_spl_cnd_rd_wnnumap_plot_png
    - umap_spl_ph_rd_rnaumap_plot_png
    - umap_spl_ph_rd_atacumap_plot_png
    - umap_spl_ph_rd_wnnumap_plot_png
    - cmp_gr_ctyp_spl_idnt_plot_png
    - cmp_gr_idnt_spl_ctyp_plot_png
    - cmp_gr_ph_spl_idnt_plot_png
    - cmp_gr_ctyp_spl_cnd_plot_png
    - cmp_gr_cnd_spl_ctyp_plot_png
    - cmp_gr_ph_spl_ctyp_plot_png
    - xpr_avg_plot_png
    - xpr_dnst_plot_png
    - xpr_per_cell_rd_rnaumap_plot_png
    - xpr_per_cell_rd_atacumap_plot_png
    - xpr_per_cell_rd_wnnumap_plot_png
    - xpr_per_cell_sgnl_rd_rnaumap_plot_png
    - xpr_per_cell_sgnl_rd_atacumap_plot_png
    - xpr_per_cell_sgnl_rd_wnnumap_plot_png
    - cvrg_plot_png
    - xpr_htmp_plot_png
    - gene_markers_tsv
    - peak_markers_tsv
    - ucsc_cb_config_data
    - ucsc_cb_html_data
    - ucsc_cb_html_file
    - seurat_data_rds
    - stdout_log
    - stderr_log

  compress_cellbrowser_config_data:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: ctype_assign/ucsc_cb_config_data
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell Manual Cell Type Assignment"
s:name: "Single-cell Manual Cell Type Assignment"
s:alternateName: "Assigns cell types for clusters based on the provided metadata file"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-ctype-assign.cwl
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
  Single-cell Manual Cell Type Assignment
  
  Assigns cell types for clusters based on the provided metadata file.