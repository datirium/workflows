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
  - "sc-rna-cluster.cwl"
  - "sc-atac-cluster.cwl"
  - "sc-wnn-cluster.cwl"
  sc_arc_sample:
  - "cellranger-arc-count.cwl"
  - "cellranger-arc-aggr.cwl"


inputs:

  alias:
    type: string
    label: "Analysis name"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Single-cell Cluster Analysis"
    doc: |
      Analysis that includes clustered
      single-cell data and was run through
      at least one of the following workflows:
      "Single-cell RNA-Seq Cluster Analysis",
      "Single-cell ATAC-Seq Cluster Analysis",
      "Single-cell WNN Cluster Analysis", -
      at any of the processing stages.
    'sd:upstreamSource': "sc_tools_sample/seurat_data_rds"
    'sd:localLabel': true

  atac_fragments_file:
    type: File?
    secondaryFiles:
    - .tbi
    label: "Cell Ranger ARC Sample (optional)"
    doc: |
      "Cell Ranger ARC Sample" for generating
      ATAC fragments coverage plots over the genes
      of interest.
    'sd:upstreamSource': "sc_arc_sample/atac_fragments_file"
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
    label: "Dimensionality reduction"
    doc: |
      Dimensionality reduction for which
      cluster names should be assigned.
  
  query_resolution:
    type: float
    label: "Clustering resolution"
    doc: |
      Clustering resolution for the selected
      "Dimensionality reduction" to be used
      for cluster names assignment.

  query_splitby_column:
    type:
    - "null"
    - type: enum
      symbols:
      - "dataset"
      - "condition"
      - "none"
    default: "none"
    label: "Criteria to split every cluster by (optional)"
    doc: |
      Criteria to split every cluster defined by
      the selected dimensionality reduction and
      resolution into several groups.
      Default: "none"

  identify_diff_genes:
    type: boolean?
    default: true
    label: "Find gene markers"
    doc: |
      Identify upregulated genes in each
      cell type compared to all other cells.
      Include only genes that are expressed
      in at least 10% of the cells coming
      from either current cell type or from
      all other cell types together.
      Exclude cells with log2FoldChange
      values less than 0.25. Use Wilcoxon
      Rank Sum test to calculate P-values.
      Keep only genes with P-values lower
      than 0.01. Adjust P-values for multiple
      comparisons using Bonferroni correction.
      Default: true

  identify_diff_peaks:
    type: boolean?
    default: false
    label: "Find peak markers"
    doc: |
      Identify differentially accessible
      peaks in each cell type compared to
      all other cells. Include only peaks
      that are present in at least 5% of
      the cells coming from either current
      cell type or from all other cell
      types together. Exclude cells with
      log2FoldChange values less than 0.25.
      Use logistic regression framework to
      calculate P-values. Keep only genes
      with P-values lower than 0.01. Adjust
      P-values for multiple comparisons
      using Bonferroni correction.
      Default: false

  genes_of_interest:
    type: string?
    default: null
    label: "Genes of interest"
    doc: |
      Comma or space separated list of genes
      of interest to visualize expression and
      to generate ATAC fragments coverage plots.
      Ignored if "Cell Ranger ARC Sample" input
      is not provided.
      Default: None

  cell_type_data:
    type: File
    label: "Cell types"
    doc: |
      A TSV/CSV file with the names for each
      cluster defined by "Clustering resolution"
      and "Dimensionality reduction" parameters.
      The file should have two columns named
      'cluster' and 'celltype'.

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
    label: "Plots color theme"
    doc: |
      Color theme for all plots saved
      as PNG files.
      Default: classic
    "sd:layout":
      advanced: true

  threads:
    type:
    - "null"
    - type: enum
      symbols:
      - "1"
      - "2"
      - "3"
      - "4"
      - "5"
      - "6"
    default: "1"
    label: "Cores/CPUs"
    doc: |
      Parallelization parameter to define the
      number of cores/CPUs that can be utilized
      simultaneously.
      Default: 1
    "sd:layout":
      advanced: true


outputs:

  umap_rd_rnaumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_rd_rnaumap_plot_png
    label: "UMAP, colored by cell type, RNA"
    doc: |
      UMAP, colored by cell type, RNA
    'sd:visualPlugins':
    - image:
        tab: 'Per cell type'
        Caption: 'UMAP, colored by cell type, RNA'

  umap_rd_atacumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_rd_atacumap_plot_png
    label: "UMAP, colored by cell type, ATAC"
    doc: |
      UMAP, colored by cell type, ATAC
    'sd:visualPlugins':
    - image:
        tab: 'Per cell type'
        Caption: 'UMAP, colored by cell type, ATAC'

  umap_rd_wnnumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_rd_wnnumap_plot_png
    label: "UMAP, colored by cell type, WNN"
    doc: |
      UMAP, colored by cell type, WNN
    'sd:visualPlugins':
    - image:
        tab: 'Per cell type'
        Caption: 'UMAP, colored by cell type, WNN'

  umap_spl_ph_rd_rnaumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_ph_rd_rnaumap_plot_png
    label: "UMAP, colored by cell type, split by cell cycle phase, RNA"
    doc: |
      UMAP, colored by cell type, split
      by cell cycle phase, RNA
    'sd:visualPlugins':
    - image:
        tab: 'Per cell type'
        Caption: 'UMAP, colored by cell type, split by cell cycle phase, RNA'

  umap_spl_ph_rd_atacumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_ph_rd_atacumap_plot_png
    label: "UMAP, colored by cell type, split by cell cycle phase, ATAC"
    doc: |
      UMAP, colored by cell type, split
      by cell cycle phase, ATAC
    'sd:visualPlugins':
    - image:
        tab: 'Per cell type'
        Caption: 'UMAP, colored by cell type, split by cell cycle phase, ATAC'

  umap_spl_ph_rd_wnnumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_ph_rd_wnnumap_plot_png
    label: "UMAP, colored by cell type, split by cell cycle phase, WNN"
    doc: |
      UMAP, colored by cell type, split
      by cell cycle phase, WNN
    'sd:visualPlugins':
    - image:
        tab: 'Per cell type'
        Caption: 'UMAP, colored by cell type, split by cell cycle phase, WNN'

  cmp_gr_ph_spl_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_ph_spl_ctyp_plot_png
    label: "Composition plot, colored by cell cycle phase, split by cell type, downsampled"
    doc: |
      Composition plot, colored by cell
      cycle phase, split by cell type,
      downsampled
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Composition plot, colored by cell cycle phase, split by cell type, downsampled'

  umap_spl_idnt_rd_rnaumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_idnt_rd_rnaumap_plot_png
    label: "UMAP, colored by cell type, split by dataset, RNA"
    doc: |
      UMAP, colored by cell type,
      split by dataset, RNA
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'UMAP, colored by cell type, split by dataset, RNA'

  umap_spl_idnt_rd_atacumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_idnt_rd_atacumap_plot_png
    label: "UMAP, colored by cell type, split by dataset, ATAC"
    doc: |
      UMAP, colored by cell type,
      split by dataset, ATAC
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'UMAP, colored by cell type, split by dataset, ATAC'

  umap_spl_idnt_rd_wnnumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_idnt_rd_wnnumap_plot_png
    label: "UMAP, colored by cell type, split by dataset, WNN"
    doc: |
      UMAP, colored by cell type,
      split by dataset, WNN
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'UMAP, colored by cell type, split by dataset, WNN'

  cmp_gr_ctyp_spl_idnt_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_ctyp_spl_idnt_plot_png
    label: "Composition plot, colored by cell type, split by dataset, downsampled"
    doc: |
      Composition plot, colored by cell
      type, split by dataset, downsampled
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Composition plot, colored by cell type, split by dataset, downsampled'

  cmp_gr_idnt_spl_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_idnt_spl_ctyp_plot_png
    label: "Composition plot, colored by dataset, split by cell type, downsampled"
    doc: |
      Composition plot, colored by
      dataset, split by cell type,
      downsampled
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Composition plot, colored by dataset, split by cell type, downsampled'

  cmp_gr_ph_spl_idnt_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_ph_spl_idnt_plot_png
    label: "Composition plot, colored by cell cycle phase, split by dataset, downsampled"
    doc: |
      Composition plot, colored by
      cell cycle phase, split by
      dataset, downsampled
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Composition plot, colored by cell cycle phase, split by dataset, downsampled'

  umap_spl_cnd_rd_rnaumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_cnd_rd_rnaumap_plot_png
    label: "UMAP, colored by cell type, split by grouping condition, RNA"
    doc: |
      UMAP, colored by cell type, split
      by grouping condition, RNA
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'UMAP, colored by cell type, split by grouping condition, RNA'

  umap_spl_cnd_rd_atacumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_cnd_rd_atacumap_plot_png
    label: "UMAP, colored by cell type, split by grouping condition, ATAC"
    doc: |
      UMAP, colored by cell type, split
      by grouping condition, ATAC
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'UMAP, colored by cell type, split by grouping condition, ATAC'

  umap_spl_cnd_rd_wnnumap_plot_png:
    type: File?
    outputSource: ctype_assign/umap_spl_cnd_rd_wnnumap_plot_png
    label: "UMAP, colored by cell type, split by grouping condition, WNN"
    doc: |
      UMAP, colored by cell type, split
      by grouping condition, WNN
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'UMAP, colored by cell type, split by grouping condition, WNN'

  cmp_gr_ctyp_spl_cnd_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_ctyp_spl_cnd_plot_png
    label: "Composition plot, colored by cell type, split by grouping condition, downsampled"
    doc: |
      Composition plot, colored by cell
      type, split by grouping condition,
      downsampled
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Composition plot, colored by cell type, split by grouping condition, downsampled'

  cmp_gr_cnd_spl_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_cnd_spl_ctyp_plot_png
    label: "Composition plot, colored by grouping condition, split by cell type, downsampled"
    doc: |
      Composition plot, colored by
      grouping condition, split by
      cell type, downsampled
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Composition plot, colored by grouping condition, split by cell type, downsampled'

  xpr_avg_plot_png:
    type: File?
    outputSource: ctype_assign/xpr_avg_plot_png
    label: "Gene expression dot plot"
    doc: |
      Gene expression dot plot
    'sd:visualPlugins':
    - image:
        tab: 'Genes of interest'
        Caption: 'Gene expression dot plot'

  xpr_dnst_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: ctype_assign/xpr_dnst_plot_png
    label: "Gene expression violin plot"
    doc: |
      Gene expression violin plot
    'sd:visualPlugins':
    - image:
        tab: 'Genes of interest'
        Caption: 'Gene expression violin plot'

  xpr_per_cell_rd_rnaumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: ctype_assign/xpr_per_cell_rd_rnaumap_plot_png
    label: "UMAP, gene expression, RNA"
    doc: |
      UMAP, gene expression, RNA
    'sd:visualPlugins':
    - image:
        tab: 'Genes of interest'
        Caption: 'UMAP, gene expression, RNA'

  xpr_per_cell_rd_atacumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: ctype_assign/xpr_per_cell_rd_atacumap_plot_png
    label: "UMAP, gene expression, ATAC"
    doc: |
      UMAP, gene expression, ATAC
    'sd:visualPlugins':
    - image:
        tab: 'Genes of interest'
        Caption: 'UMAP, gene expression, ATAC'

  xpr_per_cell_rd_wnnumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: ctype_assign/xpr_per_cell_rd_wnnumap_plot_png
    label: "UMAP, gene expression, WNN"
    doc: |
      UMAP, gene expression, WNN
    'sd:visualPlugins':
    - image:
        tab: 'Genes of interest'
        Caption: 'UMAP, gene expression, WNN'

  xpr_htmp_plot_png:
    type: File?
    outputSource: ctype_assign/xpr_htmp_plot_png
    label: "Gene expression heatmap"
    doc: |
      Gene expression heatmap
    'sd:visualPlugins':
    - image:
        tab: 'Heatmap'
        Caption: 'Gene expression heatmap'

  cvrg_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: ctype_assign/cvrg_plot_png
    label: "ATAC fragments coverage"
    doc: |
      ATAC fragments coverage
    'sd:visualPlugins':
    - image:
        tab: 'Genome coverage'
        Caption: 'ATAC fragments coverage'

  xpr_htmp_tsv:
    type: File?
    outputSource: ctype_assign/xpr_htmp_tsv
    label: "Markers from gene expression heatmap"
    doc: |
      Gene markers used for gene
      expression heatmap

  gene_markers_tsv:
    type: File?
    outputSource: ctype_assign/gene_markers_tsv
    label: "Gene markers per cell type"
    doc: |
      Gene markers per cell type
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Gene markers'
        Title: 'Gene markers per cell type'

  peak_markers_tsv:
    type: File?
    outputSource: ctype_assign/peak_markers_tsv
    label: "Peak markers per cell type"
    doc: |
      Peak markers per cell type
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Peak markers'
        Title: 'Peak markers per cell type'

  ucsc_cb_html_data:
    type: Directory?
    outputSource: ctype_assign/ucsc_cb_html_data
    label: "UCSC Cell Browser data"
    doc: |
      Directory with UCSC Cell Browser
      data

  ucsc_cb_html_file:
    type: File?
    outputSource: ctype_assign/ucsc_cb_html_file
    label: "UCSC Cell Browser"
    doc: |
      UCSC Cell Browser HTML index file
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  seurat_data_rds:
    type: File
    outputSource: ctype_assign/seurat_data_rds
    label: "Processed Seurat data in RDS format"
    doc: |
      Processed Seurat data in RDS format

  seurat_data_scope:
    type: File?
    outputSource: ctype_assign/seurat_data_scope
    label: "Processed Seurat data in SCope compatible loom format"
    doc: |
      Processed Seurat data in SCope compatible loom format

  pdf_plots:
    type: File
    outputSource: compress_pdf_plots/compressed_folder
    label: "Plots in PDF format"
    doc: |
      Compressed folder with plots
      in PDF format

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
      query_splitby_column:
        source: query_splitby_column
        valueFrom: |
          ${
            if (self == "dataset") {
              return "new.ident";
            } else if (self == "condition") {
              return "condition";
            } else {
              return null;
            }
          }
      atac_fragments_file: atac_fragments_file
      genes_of_interest:
        source: genes_of_interest
        valueFrom: $(split_features(self))
      identify_diff_genes: identify_diff_genes
      identify_diff_peaks: identify_diff_peaks
      rna_minimum_logfc:
        default: 0.25
      rna_minimum_pct:
        default: 0.1
      atac_minimum_logfc:
        default: 0.25
      atac_minimum_pct:
        default: 0.05
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
      export_scope_data:
        default: true
      export_pdf_plots:
        default: true
      color_theme: color_theme
      parallel_memory_limit:
        default: 32
      vector_memory_limit:
        default: 96
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
    - cvrg_plot_png
    - xpr_htmp_plot_png
    - umap_rd_rnaumap_plot_pdf
    - umap_rd_atacumap_plot_pdf
    - umap_rd_wnnumap_plot_pdf
    - umap_spl_idnt_rd_rnaumap_plot_pdf
    - umap_spl_idnt_rd_atacumap_plot_pdf
    - umap_spl_idnt_rd_wnnumap_plot_pdf
    - umap_spl_cnd_rd_rnaumap_plot_pdf
    - umap_spl_cnd_rd_atacumap_plot_pdf
    - umap_spl_cnd_rd_wnnumap_plot_pdf
    - umap_spl_ph_rd_rnaumap_plot_pdf
    - umap_spl_ph_rd_atacumap_plot_pdf
    - umap_spl_ph_rd_wnnumap_plot_pdf
    - cmp_gr_ctyp_spl_idnt_plot_pdf
    - cmp_gr_idnt_spl_ctyp_plot_pdf
    - cmp_gr_ph_spl_idnt_plot_pdf
    - cmp_gr_ctyp_spl_cnd_plot_pdf
    - cmp_gr_cnd_spl_ctyp_plot_pdf
    - cmp_gr_ph_spl_ctyp_plot_pdf
    - xpr_avg_plot_pdf
    - xpr_dnst_plot_pdf
    - xpr_per_cell_rd_rnaumap_plot_pdf
    - xpr_per_cell_rd_atacumap_plot_pdf
    - xpr_per_cell_rd_wnnumap_plot_pdf
    - xpr_per_cell_sgnl_rd_rnaumap_plot_pdf
    - xpr_per_cell_sgnl_rd_atacumap_plot_pdf
    - xpr_per_cell_sgnl_rd_wnnumap_plot_pdf
    - cvrg_plot_pdf
    - xpr_htmp_plot_pdf
    - xpr_htmp_tsv
    - gene_markers_tsv
    - peak_markers_tsv
    - ucsc_cb_html_data
    - ucsc_cb_html_file
    - seurat_data_rds
    - seurat_data_scope
    - stdout_log
    - stderr_log

  pdf_plots:
    run: ../tools/files-to-folder.cwl
    in:
      input_files:
        source:
        - ctype_assign/umap_rd_rnaumap_plot_pdf
        - ctype_assign/umap_rd_atacumap_plot_pdf
        - ctype_assign/umap_rd_wnnumap_plot_pdf
        - ctype_assign/umap_spl_idnt_rd_rnaumap_plot_pdf
        - ctype_assign/umap_spl_idnt_rd_atacumap_plot_pdf
        - ctype_assign/umap_spl_idnt_rd_wnnumap_plot_pdf
        - ctype_assign/umap_spl_cnd_rd_rnaumap_plot_pdf
        - ctype_assign/umap_spl_cnd_rd_atacumap_plot_pdf
        - ctype_assign/umap_spl_cnd_rd_wnnumap_plot_pdf
        - ctype_assign/umap_spl_ph_rd_rnaumap_plot_pdf
        - ctype_assign/umap_spl_ph_rd_atacumap_plot_pdf
        - ctype_assign/umap_spl_ph_rd_wnnumap_plot_pdf
        - ctype_assign/cmp_gr_ctyp_spl_idnt_plot_pdf
        - ctype_assign/cmp_gr_idnt_spl_ctyp_plot_pdf
        - ctype_assign/cmp_gr_ph_spl_idnt_plot_pdf
        - ctype_assign/cmp_gr_ctyp_spl_cnd_plot_pdf
        - ctype_assign/cmp_gr_cnd_spl_ctyp_plot_pdf
        - ctype_assign/cmp_gr_ph_spl_ctyp_plot_pdf
        - ctype_assign/xpr_avg_plot_pdf
        - ctype_assign/xpr_dnst_plot_pdf
        - ctype_assign/xpr_per_cell_rd_rnaumap_plot_pdf
        - ctype_assign/xpr_per_cell_rd_atacumap_plot_pdf
        - ctype_assign/xpr_per_cell_rd_wnnumap_plot_pdf
        - ctype_assign/xpr_per_cell_sgnl_rd_rnaumap_plot_pdf
        - ctype_assign/xpr_per_cell_sgnl_rd_atacumap_plot_pdf
        - ctype_assign/xpr_per_cell_sgnl_rd_wnnumap_plot_pdf
        - ctype_assign/cvrg_plot_pdf
        - ctype_assign/xpr_htmp_plot_pdf
        valueFrom: $(self.flat().filter(n => n))
      folder_basename:
        default: "pdf_plots"
    out:
    - folder

  compress_pdf_plots:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: pdf_plots/folder
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

  Assigns cell types for clusters based on
  the provided metadata file.