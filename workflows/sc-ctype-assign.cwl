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


"sd:upstream":
  sc_tools_sample:
  - "sc-rna-cluster.cwl"
  - "sc-atac-cluster.cwl"
  - "sc-wnn-cluster.cwl"
  - "sc-rna-azimuth.cwl"
  sc_atac_sample:
  - "cellranger-arc-count.cwl"
  - "cellranger-arc-aggr.cwl"
  - "cellranger-atac-count.cwl"
  - "cellranger-atac-aggr.cwl"


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
      "Single-Cell RNA-Seq Cluster Analysis",
      "Single-Cell ATAC-Seq Cluster Analysis",
      "Single-Cell WNN Cluster Analysis", -
      at any of the processing stages.
    "sd:upstreamSource": "sc_tools_sample/seurat_data_rds"
    "sd:localLabel": true

  atac_fragments_file:
    type: File?
    secondaryFiles:
    - .tbi
    label: "Cell Ranger ATAC or RNA+ATAC Sample (optional)"
    doc: |
      Any "Cell Ranger ATAC or RNA+ATAC Sample"
      for generating ATAC fragments coverage
      plots over the genes of interest. This
      sample can be obtained from one of the
      following pipelines: "Cell Ranger Count
      (RNA+ATAC)", "Cell Ranger Aggregate
      (RNA+ATAC)", "Cell Ranger Count (ATAC)",
      or "Cell Ranger Aggregate (ATAC)".
    "sd:upstreamSource": "sc_atac_sample/atac_fragments_file"
    "sd:localLabel": true

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
      Ignored if "Cell Ranger ATAC or RNA+ATAC
      Sample (optional)" input is not provided.
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

  barcodes_data:
    type: File?
    label: "Selected cell barcodes (optional)"
    doc: |
      A TSV/CSV file to optionally extend
      the single cell metadata with the custom
      values per each barcode. The provided
      file should include the first column
      named "barcode", with one cell barcode
      per line. All other columns, except for
      the "barcode", will be added to the single
      cell metadata loaded from the "Single-cell
      Cluster Analysis" and can be utilized in
      the current or future steps of analysis.

  genesets_data:
    type: File?
    label: "GMT file for calculating average expression levels per gene set (optional)"
    doc: |
      Path to the GMT file for calculating average expression levels
      (module scores) per gene set. This file can be downloaded from
      the Molecular Signatures Database (MSigDB) following the link
      https://www.gsea-msigdb.org/gsea/msigdb. To calculate module
      scores the loaded Seurat object should include RNA assay.
      Default: do not calculate gene set expression scores.

  export_loupe_data:
    type: boolean?
    default: false
    label: "Save raw counts to Loupe file. I confirm that data is generated by 10x technology and accept the EULA available at https://10xgen.com/EULA"
    doc: |
      Save raw counts from the RNA assay to Loupe file. By
      enabling this feature you accept the End-User License
      Agreement available at https://10xgen.com/EULA.
      Default: false
    "sd:layout":
      advanced: true

  export_html_report:
    type: boolean?
    default: true
    label: "Show HTML report"
    doc: |
      Export tehcnical report in HTML format.
      Default: true
    "sd:layout":
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
    default: "4"
    label: "Cores/CPUs"
    doc: |
      Parallelization parameter to define the
      number of cores/CPUs that can be utilized
      simultaneously.
      Default: 4
    "sd:layout":
      advanced: true


outputs:

  cell_cnts_gr_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/cell_cnts_gr_ctyp_plot_png
    label: "Number of cells per cell type (all cells)"
    doc: |
      Number of cells per cell type.
      All cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Number of cells per cell type (all cells)"

  qc_mtrcs_dnst_gr_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/qc_mtrcs_dnst_gr_ctyp_plot_png
    label: "Distribution of QC metrics per cell colored by cell type (all cells)"
    doc: |
      Distribution of QC metrics per cell
      colored by cell type.
      All cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Distribution of QC metrics per cell colored by cell type (all cells)"

  gene_umi_spl_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/gene_umi_spl_ctyp_plot_png
    label: "Genes vs RNA reads per cell (split by cell type, all cells)"
    doc: |
      Genes vs RNA reads per cell.
      Split by cell type; all cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Genes vs RNA reads per cell (split by cell type, all cells)"

  umi_mito_spl_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/umi_mito_spl_ctyp_plot_png
    label: "RNA reads vs mitochondrial % per cell (split by cell type, all cells)"
    doc: |
      RNA reads vs mitochondrial % per cell.
      Split by cell type; all cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "RNA reads vs mitochondrial % per cell (split by cell type, all cells)"

  tss_frgm_spl_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/tss_frgm_spl_ctyp_plot_png
    label: "TSS enrichment score vs ATAC fragments in peaks per cell (split by cell type, all cells)"
    doc: |
      TSS enrichment score vs ATAC
      fragments in peaks per cell.
      Split by cell type; all cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "TSS enrichment score vs ATAC fragments in peaks per cell (split by cell type, all cells)"

  rna_atac_cnts_spl_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/rna_atac_cnts_spl_ctyp_plot_png
    label: "RNA reads vs ATAC fragments in peaks per cell (split by cell type, all cells)"
    doc: |
      RNA reads vs ATAC fragments in peaks per cell.
      Split by cell type; all cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "RNA reads vs ATAC fragments in peaks per cell (split by cell type, all cells)"

  rnadbl_gr_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/rnadbl_gr_ctyp_plot_png
    label: "Percentage of RNA doublets per cell type (all cells)"
    doc: |
      Percentage of RNA doublets per cell type.
      All cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Percentage of RNA doublets per cell type (all cells)"

  atacdbl_gr_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/atacdbl_gr_ctyp_plot_png
    label: "Percentage of ATAC doublets per cell type (all cells)"
    doc: |
      Percentage of ATAC doublets per cell type.
      All cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Percentage of ATAC doublets per cell type (all cells)"

  vrlpdbl_gr_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/vrlpdbl_gr_ctyp_plot_png
    label: "Percentage of RNA and ATAC doublets per cell type (all cells)"
    doc: |
      Percentage of RNA and ATAC doublets
      per cell type.
      All cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Percentage of RNA and ATAC doublets per cell type (all cells)"

  umap_gr_ctyp_plot_png:
    type: File?
    outputSource: ctype_assign/umap_gr_ctyp_plot_png
    label: "UMAP colored by cell type (all cells)"
    doc: |
      UMAP colored by cell type.
      All cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by cell type"
        Caption: "UMAP colored by cell type (all cells)"

  umap_gr_ctyp_spl_ph_png:
    type: File?
    outputSource: ctype_assign/umap_gr_ctyp_spl_ph_png
    label: "UMAP colored by cell type (split by cell cycle phase, optionally downsampled)"
    doc: |
      UMAP colored by cell type.
      Split by cell cycle phase; downsampled
      to the smallest dataset (if multiple
      datasets are analyzed jointly).
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by cell type"
        Caption: "UMAP colored by cell type (split by cell cycle phase, optionally downsampled)"

  cmp_gr_ph_spl_ctyp_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_ph_spl_ctyp_png
    label: "Composition plot colored by cell cycle phase (split by cell type, optionally downsampled)"
    doc: |
      Composition plot colored by cell cycle phase.
      Split by cell type; downsampled to the
      smallest dataset (if multiple datasets are
      analyzed jointly).
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by cell type"
        Caption: "Composition plot colored by cell cycle phase (split by cell type, optionally downsampled)"

  umap_gr_ctyp_spl_idnt_plot_png:
    type: File?
    outputSource: ctype_assign/umap_gr_ctyp_spl_idnt_plot_png
    label: "UMAP colored by cell type (split by dataset, downsampled)"
    doc: |
      UMAP colored by cell type.
      Split by dataset; downsampled to the
      smallest dataset.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by dataset"
        Caption: "UMAP colored by cell type (split by dataset, downsampled)"

  cmp_gr_ctyp_spl_idnt_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_ctyp_spl_idnt_plot_png
    label: "Composition plot colored by cell type (split by dataset, downsampled)"
    doc: |
      Composition plot colored by cell type.
      Split by dataset; downsampled to the
      smallest dataset.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by dataset"
        Caption: "Composition plot colored by cell type (split by dataset, downsampled)"

  umap_gr_ph_spl_idnt_plot_png:
    type: File?
    outputSource: ctype_assign/umap_gr_ph_spl_idnt_plot_png
    label: "UMAP colored by cell cycle phase (split by dataset, downsampled)"
    doc: |
      UMAP colored by cell cycle phase.
      Split by dataset; downsampled to the
      smallest dataset.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by dataset"
        Caption: "UMAP colored by cell cycle phase (split by dataset, downsampled)"

  cmp_gr_ph_spl_idnt_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_ph_spl_idnt_plot_png
    label: "Composition plot colored by cell cycle phase (split by dataset, downsampled)"
    doc: |
      Composition plot colored by cell cycle phase.
      Split by dataset; downsampled to the smallest
      dataset.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by dataset"
        Caption: "Composition plot colored by cell cycle phase (split by dataset, downsampled)"

  umap_gr_ctyp_spl_cnd_plot_png:
    type: File?
    outputSource: ctype_assign/umap_gr_ctyp_spl_cnd_plot_png
    label: "UMAP colored by cell type (split by grouping condition, downsampled)"
    doc: |
      UMAP colored by cell type.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by group"
        Caption: "UMAP colored by cell type (split by grouping condition, downsampled)"

  cmp_gr_ctyp_spl_cnd_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_ctyp_spl_cnd_plot_png
    label: "Composition plot colored by cell type (split by grouping condition, downsampled)"
    doc: |
      Composition plot colored by cell type.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by group"
        Caption: "Composition plot colored by cell type (split by grouping condition, downsampled)"

  umap_gr_ph_spl_cnd_plot_png:
    type: File?
    outputSource: ctype_assign/umap_gr_ph_spl_cnd_plot_png
    label: "UMAP colored by cell cycle phase (split by grouping condition, downsampled)"
    doc: |
      UMAP colored by cell cycle phase.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by group"
        Caption: "UMAP colored by cell cycle phase (split by grouping condition, downsampled)"

  cmp_gr_ph_spl_cnd_plot_png:
    type: File?
    outputSource: ctype_assign/cmp_gr_ph_spl_cnd_plot_png
    label: "Composition plot colored by cell cycle phase (split by grouping condition, downsampled)"
    doc: |
      Composition plot colored by cell cycle phase.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by group"
        Caption: "Composition plot colored by cell cycle phase (split by grouping condition, downsampled)"

  xpr_avg_plot_png:
    type: File?
    outputSource: ctype_assign/xpr_avg_plot_png
    label: "Average gene expression"
    doc: |
      Average gene expression.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Genes of interest (expression)"
        Caption: "Average gene expression"

  xpr_dnst_plot_png:
    type: File?
    outputSource: ctype_assign/xpr_dnst_plot_png
    label: "Gene expression density"
    doc: |
      Gene expression density.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Genes of interest (expression)"
        Caption: "Gene expression density"

  xpr_per_cell_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: ctype_assign/xpr_per_cell_plot_png
    label: "UMAP colored by gene expression (per gene)"
    doc: |
      UMAP colored by gene expression.
      All genes of interest.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Genes of interest (expression)"
        Caption: "UMAP colored by gene expression (per gene)"

  cvrg_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: ctype_assign/cvrg_plot_png
    label: "ATAC fragment coverage (per gene)"
    doc: |
      ATAC fragment coverage.
      All genes of interest.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Genes of interest (coverage)"
        Caption: "ATAC fragment coverage (per gene)"

  gse_per_cell_plot_png:
    type: File?
    outputSource: ctype_assign/gse_per_cell_plot_png
    label: "UMAP colored by gene set expression score"
    doc: |
      UMAP colored by gene set expression score.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Gene sets of interest (expression)"
        Caption: "UMAP colored by gene set expression score"

  gse_avg_plot_png:
    type: File?
    outputSource: ctype_assign/gse_avg_plot_png
    label: "Average gene set expression score"
    doc: |
      Average gene set expression score.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Gene sets of interest (expression)"
        Caption: "Average gene set expression score"

  gse_dnst_plot_png:
    type: File?
    outputSource: ctype_assign/gse_dnst_plot_png
    label: "Gene set expression score density"
    doc: |
      Gene set expression score density.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Gene sets of interest (expression)"
        Caption: "Gene set expression score density"

  xpr_htmp_plot_png:
    type: File?
    outputSource: ctype_assign/xpr_htmp_plot_png
    label: "Gene expression heatmap (top gene markers)"
    doc: |
      Gene expression heatmap.
      Top gene markers.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Gene markers heatmap"
        Caption: "Gene expression heatmap (top gene markers)"

  xpr_htmp_tsv:
    type: File?
    outputSource: ctype_assign/xpr_htmp_tsv
    label: "Gene expression heatmap (top gene markers)"
    doc: |
      Gene expression heatmap.
      Top gene markers.
      TSV format.

  gene_markers_tsv:
    type: File?
    outputSource: ctype_assign/gene_markers_tsv
    label: "Gene markers"
    doc: |
      Gene markers.
      TSV format.
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Gene markers table"
        Title: "Gene markers"

  peak_markers_tsv:
    type: File?
    outputSource: ctype_assign/peak_markers_tsv
    label: "Peak markers"
    doc: |
      Peak markers.
      TSV format.
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Peak markers table"
        Title: "Peak markers"

  ucsc_cb_html_data:
    type: Directory?
    outputSource: ctype_assign/ucsc_cb_html_data
    label: "UCSC Cell Browser (data)"
    doc: |
      UCSC Cell Browser html data.

  ucsc_cb_html_file:
    type: File?
    outputSource: ctype_assign/ucsc_cb_html_file
    label: "UCSC Cell Browser"
    doc: |
      UCSC Cell Browser html index.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  seurat_data_rds:
    type: File
    outputSource: ctype_assign/seurat_data_rds
    label: "Seurat object in RDS format"
    doc: |
      Seurat object.
      RDS format.

  seurat_data_scope:
    type: File?
    outputSource: ctype_assign/seurat_data_scope
    label: "Seurat object in SCope compatible loom format"
    doc: |
      Seurat object.
      SCope compatible.
      Loom format.

  reference_data_rds:
    type: File?
    outputSource: ctype_assign/reference_data_rds
    label: "Seurat object formatted as an Azimuth reference model"
    doc: |
      Seurat object with assigned cell
      types formatted as an Azimuth
      reference model.
      RDS format.

  reference_data_index:
    type: File?
    outputSource: ctype_assign/reference_data_index
    label: "Annoy index for the Azimuth reference model"
    doc: |
      Annoy index generated for the
      Azimuth reference model.
      Annoy format.

  seurat_rna_data_cloupe:
    type: File?
    outputSource: ctype_assign/seurat_rna_data_cloupe
    label: "Seurat object in Loupe format"
    doc: |
      Seurat object.
      RNA counts.
      Loupe format.

  pdf_plots:
    type: File
    outputSource: compress_pdf_plots/compressed_folder
    label: "Compressed folder with all PDF plots"
    doc: |
      Compressed folder with all PDF plots.

  sc_report_html_file:
    type: File?
    outputSource: ctype_assign/sc_report_html_file
    label: "Analysis log"
    doc: |
      Tehcnical report.
      HTML format.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  ctype_assign_human_log:
    type: File
    outputSource: ctype_assign/human_log
    label: "Human readable error log"
    doc: |
      Human readable error log
      from the ctype_assign step.

  ctype_assign_stdout_log:
    type: File
    outputSource: ctype_assign/stdout_log
    label: "Output log"
    doc: |
      Stdout log from the ctype_assign step.

  ctype_assign_stderr_log:
    type: File
    outputSource: ctype_assign/stderr_log
    label: "Error log"
    doc: |
      Stderr log from the ctype_assign step.


steps:

  ctype_assign:
    run: ../tools/sc-ctype-assign.cwl
    in:
      query_data_rds: query_data_rds
      cell_type_data: cell_type_data
      barcodes_data: barcodes_data
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
      reduction:
        source: query_reduction
        valueFrom: |
          ${
            if (self == "RNA") {
              return "rnaumap";
            } else if (self == "ATAC") {
              return "atacumap";
            } else {
              return "wnnumap";
            }
          }
      atac_fragments_file: atac_fragments_file
      genes_of_interest:
        source: genes_of_interest
        valueFrom: $(split_features(self))
      genesets_data: genesets_data
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
      export_azimuth_ref:
        default: true
      export_ucsc_cb:
        default: true
      export_scope_data:
        default: true
      export_loupe_data: export_loupe_data
      export_pdf_plots:
        default: true
      color_theme: color_theme
      parallel_memory_limit:
        default: 32
      vector_memory_limit:
        default: 128
      export_html_report: export_html_report
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - cell_cnts_gr_ctyp_plot_png
    - gene_umi_spl_ctyp_plot_png
    - umi_mito_spl_ctyp_plot_png
    - rnadbl_gr_ctyp_plot_png
    - tss_frgm_spl_ctyp_plot_png
    - atacdbl_gr_ctyp_plot_png
    - rna_atac_cnts_spl_ctyp_plot_png
    - vrlpdbl_gr_ctyp_plot_png
    - qc_mtrcs_dnst_gr_ctyp_plot_png
    - umap_gr_ctyp_plot_png
    - umap_gr_ctyp_spl_idnt_plot_png
    - cmp_gr_ctyp_spl_idnt_plot_png
    - umap_gr_ph_spl_idnt_plot_png
    - cmp_gr_ph_spl_idnt_plot_png
    - umap_gr_ctyp_spl_ph_png
    - cmp_gr_ph_spl_ctyp_png
    - umap_gr_ctyp_spl_cnd_plot_png
    - cmp_gr_ctyp_spl_cnd_plot_png
    - umap_gr_ph_spl_cnd_plot_png
    - cmp_gr_ph_spl_cnd_plot_png
    - gse_per_cell_plot_png
    - gse_avg_plot_png
    - gse_dnst_plot_png
    - xpr_avg_plot_png
    - xpr_per_cell_plot_png
    - xpr_dnst_plot_png
    - xpr_htmp_plot_png
    - cvrg_plot_png
    - all_plots_pdf
    - xpr_htmp_tsv
    - gene_markers_tsv
    - peak_markers_tsv
    - ucsc_cb_html_data
    - ucsc_cb_html_file
    - seurat_data_rds
    - seurat_data_scope
    - reference_data_rds
    - reference_data_index
    - seurat_rna_data_cloupe
    - sc_report_html_file
    - human_log
    - stdout_log
    - stderr_log

  folder_pdf_plots:
    run: ../tools/files-to-folder.cwl
    in:
      input_files:
        source:
        - ctype_assign/all_plots_pdf
        valueFrom: $(self.flat().filter(n => n))
      folder_basename:
        default: "pdf_plots"
    out:
    - folder

  compress_pdf_plots:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: folder_pdf_plots/folder
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-Cell Manual Cell Type Assignment"
s:name: "Single-Cell Manual Cell Type Assignment"
s:alternateName: "Single-Cell Manual Cell Type Assignment"

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
  Single-Cell Manual Cell Type Assignment

  Assigns identities to cells clustered with any of the “Single-Cell
  Cluster Analysis” pipelines. For “Single-Cell RNA-Seq Cluster
  Analysis” the results of this workflow are used in the “Single-Cell
  RNA-Seq Differential Expression Analysis”, “Single-Cell RNA-Seq
  Trajectory Analysis”, and — when combined with outputs from the
  “Cell Ranger Count (RNA+VDJ)” or “Cell Ranger Aggregate (RNA, RNA+VDJ)”
  workflow — in the “Single-Cell Immune Profiling Analysis” pipeline.
  For “Single-Cell ATAC-Seq Cluster Analysis”, the results of this
  workflow are used in the “Single-Cell ATAC-Seq Differential
  Accessibility Analysis” and “Single-Cell ATAC-Seq Genome Coverage”
  pipelines. For “Single-Cell WNN Cluster Analysis”, the results of
  this workflow are used in all of the above, except the “Single-Cell
  Immune Profiling Analysis” pipeline.