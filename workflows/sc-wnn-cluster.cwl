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
          let splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };
    - var split_numbers = function(line) {
          let splitted_line = line?line.split(/[\s,]+/).map(parseFloat):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };
    - var parse_range = function(line) {
          if (line.includes("-")) {
              const parts = line.split("-");
              const start = parseFloat(parts[0].trim());
              let end, step;
              if (parts[1].includes(":")) {
                  [end, step] = parts[1].split(":").map(Number);
              } else {
                  end = parseFloat(parts[1].trim());
                  step = 0.1;
              }
              const result = [];
              for (let i = start; i <= end; i = parseFloat((i + step).toFixed(10))) {
                  result.push(parseFloat(i.toFixed(10)));
              }
              return result;
          } else {
              return [parseFloat(line)];
          }
      };


"sd:upstream":
  sc_tools_sample:
  - "sc-wnn-cluster.cwl"
  - "sc-rna-cluster.cwl"
  - "sc-atac-cluster.cwl"
  - "sc-rna-reduce.cwl"
  - "sc-atac-reduce.cwl"
  - "sc-rna-azimuth.cwl"
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
    label: "Single-cell Analysis with both PCA and LSI Transformed Datasets"
    doc: |
      Analysis that includes single-cell
      multiome RNA and ATAC-Seq datasets
      run through both "Single-Cell
      RNA-Seq Dimensionality Reduction
      Analysis" and "Single-Cell ATAC-Seq
      Dimensionality Reduction Analysis"
      at any of the processing stages.
    "sd:upstreamSource": "sc_tools_sample/seurat_data_rds"
    "sd:localLabel": true

  atac_fragments_file:
    type: File?
    secondaryFiles:
    - .tbi
    label: "Cell Ranger RNA+ATAC Sample (optional)"
    doc: |
      Any "Cell Ranger ATAC or RNA+ATAC Sample"
      for generating ATAC fragments coverage
      plots over the genes of interest. This
      sample can be obtained from either
      "Cell Ranger Count (RNA+ATAC)" or "Cell
      Ranger Aggregate (RNA+ATAC)" pipeline
    "sd:upstreamSource": "sc_arc_sample/atac_fragments_file"
    "sd:localLabel": true

  rna_dimensions:
    type: int?
    default: 40
    label: "Target RNA dimensionality"
    doc: |
      Target RNA dimensionality is the number of
      principal components to be used in
      constructing the weighted nearest-neighbor
      graph before clustering. The accepted values
      range from 1 to 50.
      Default: 40

  atac_dimensions:
    type: int?
    default: 40
    label: "Target ATAC dimensionality"
    doc: |
      Target ATAC dimensionality is the number of
      LSI dimensions to be used in constructing
      the weighted nearest-neighbor graph before
      clustering. The accepted values range from
      2 to 50. The first dimension is always
      excluded.
      Default: 40

  resolution:
    type: string?
    default: "0.3"
    label: "Clustering resolution"
    doc: |
      The resolution defines the “granularity”
      of the clustered data. Larger resolution
      values lead to more clusters. The optimal
      resolution often increases with the number
      of cells. To run the analysis with multiple
      resolutions, provide a range in a form of
      start-end:step. Step parameter is optional
      and equal to 0.1 by default.
      Default: 0.3

  identify_diff_genes:
    type: boolean?
    default: true
    label: "Find gene markers"
    doc: |
      The user can identify upregulated genes
      in each cluster compared to all other
      cells. The results include only genes
      that are expressed in at least 10% of
      the cells coming from either the current
      cluster or from all other clusters together.
      Genes with the log2FoldChange values smaller
      than 0.25 are excluded. The p-values are
      calculated with the Wilcoxon Rank Sum test
      and adjusted for multiple comparisons using
      the Bonferroni correction.
      Default: true

  identify_diff_peaks:
    type: boolean?
    default: false
    label: "Find peak markers"
    doc: |
      The user can identify differentially accessible
      peaks in each cluster compared to all other cells.
      The results include only peaks that are present
      in at least 5% of the cells coming from either
      the current cluster or from all other clusters
      together. Peaks with log2FoldChange values smaller
      than 0.25 are excluded. The p-values are calculated
      using the logistic regression framework and adjusted
      for multiple comparisons using the Bonferroni
      correction.
      Default: false

  genes_of_interest:
    type: string?
    default: null
    label: "Genes of interest"
    doc: |
      A comma- or space-separated list of genes
      of interest to visualize expression. If the
      “Cell Ranger RNA+ATAC Sample (optional)”
      input was provided the ATAC fragment coverage
      plots will be created as well.
      Default: None

  genesets_data:
    type: File?
    label: "GMT file for calculating average expression levels per gene set (optional)"
    doc: |
      Path to the GMT file for calculating average expression levels
      (module scores) per gene set. This file can be downloaded from
      the Molecular Signatures Database (MSigDB) following the link
      https://www.gsea-msigdb.org/gsea/msigdb.
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

  cell_cnts_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/cell_cnts_gr_clst_res_plot_png
    label: "Number of cells per cluster (all cells)"
    doc: |
      Number of cells per cluster.
      All cells; all resolutions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Number of cells per cluster (all cells)"

  qc_mtrcs_dnst_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/qc_mtrcs_dnst_gr_clst_res_plot_png
    label: "Distribution of QC metrics per cell colored by cluster (all cells)"
    doc: |
      Distribution of QC metrics per cell
      colored by cluster.
      All cells; all resolutions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Distribution of QC metrics per cell colored by cluster (all cells)"

  gene_umi_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/gene_umi_spl_clst_res_plot_png
    label: "Genes vs RNA reads per cell (split by cluster, all cells)"
    doc: |
      Genes vs RNA reads per cell.
      Split by cluster; all cells;
      all resolutions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Genes vs RNA reads per cell (split by cluster, all cells)"

  umi_mito_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/umi_mito_spl_clst_res_plot_png
    label: "RNA reads vs mitochondrial % per cell (split by cluster, all cells)"
    doc: |
      RNA reads vs mitochondrial % per cell.
      Split by cluster; all cells; all
      resolutions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "RNA reads vs mitochondrial % per cell (split by cluster, all cells)"

  tss_frgm_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/tss_frgm_spl_clst_res_plot_png
    label: "TSS enrichment score vs ATAC fragments in peaks per cell (split by cluster, all cells)"
    doc: |
      TSS enrichment score vs ATAC
      fragments in peaks per cell.
      Split by cluster; all cells;
      all resolutions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "TSS enrichment score vs ATAC fragments in peaks per cell (split by cluster, all cells)"

  rna_atac_cnts_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/rna_atac_cnts_spl_clst_res_plot_png
    label: "RNA reads vs ATAC fragments in peaks per cell (split by cluster, all cells)"
    doc: |
      RNA reads vs ATAC fragments in peaks per cell.
      Split by cluster; all cells; all resolutions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "RNA reads vs ATAC fragments in peaks per cell (split by cluster, all cells)"

  rnadbl_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/rnadbl_gr_clst_res_plot_png
    label: "Percentage of RNA doublets per cluster (all cells)"
    doc: |
      Percentage of RNA doublets per cluster.
      All cells; all resolutions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Percentage of RNA doublets per cluster (all cells)"

  atacdbl_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/atacdbl_gr_clst_res_plot_png
    label: "Percentage of ATAC doublets per cluster (all cells)"
    doc: |
      Percentage of ATAC doublets per cluster.
      All cells; all resolutions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Percentage of ATAC doublets per cluster (all cells)"

  vrlpdbl_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/vrlpdbl_gr_clst_res_plot_png
    label: "Percentage of RNA and ATAC doublets per cluster (all cells)"
    doc: |
      Percentage of RNA and ATAC doublets
      per cluster.
      All cells; all resolutions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Percentage of RNA and ATAC doublets per cluster (all cells)"

  umap_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/umap_gr_clst_res_plot_png
    label: "UMAP colored by cluster (all cells)"
    doc: |
      UMAP colored by cluster.
      All cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by cluster"
        Caption: "UMAP colored by cluster (all cells)"

  slh_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/slh_gr_clst_res_plot_png
    label: "Silhouette scores (all cells)"
    doc: |
      Silhouette scores.
      All cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by cluster"
        Caption: "Silhouette scores (all cells)"

  umap_gr_clst_spl_ph_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/umap_gr_clst_spl_ph_res_plot_png
    label: "UMAP colored by cluster (split by cell cycle phase, optionally downsampled)"
    doc: |
      UMAP colored by cluster.
      Split by cell cycle phase; downsampled
      to the smallest dataset (if multiple
      datasets are analyzed jointly).
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by cluster"
        Caption: "UMAP colored by cluster (split by cell cycle phase, optionally downsampled)"

  cmp_gr_ph_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/cmp_gr_ph_spl_clst_res_plot_png
    label: "Composition plot colored by cell cycle phase (split by cluster, optionally downsampled)"
    doc: |
      Composition plot colored by cell cycle phase.
      Split by cluster; downsampled to the smallest
      dataset (if multiple datasets are analyzed
      jointly).
      PNG format
    "sd:visualPlugins":
    - image:
        tab: "Split by cluster"
        Caption: "Composition plot colored by cell cycle phase (split by cluster, optionally downsampled)"

  umap_gr_clst_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/umap_gr_clst_spl_idnt_res_plot_png
    label: "UMAP colored by cluster (split by dataset, downsampled)"
    doc: |
      UMAP colored by cluster.
      Split by dataset; downsampled
      to the smallest dataset.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by dataset"
        Caption: "UMAP colored by cluster (split by dataset, downsampled)"

  cmp_gr_clst_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/cmp_gr_clst_spl_idnt_res_plot_png
    label: "Composition plot colored by cluster (split by dataset, downsampled)"
    doc: |
      Composition plot colored by cluster.
      Split by dataset; downsampled
      to the smallest dataset.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by dataset"
        Caption: "Composition plot colored by cluster (split by dataset, downsampled)"

  umap_gr_ph_spl_idnt_plot_png:
    type: File?
    outputSource: sc_wnn_cluster/umap_gr_ph_spl_idnt_plot_png
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
    outputSource: sc_wnn_cluster/cmp_gr_ph_spl_idnt_plot_png
    label: "Composition plot colored by cell cycle phase (split by dataset, downsampled)"
    doc: |
      Composition plot colored by cell cycle phase.
      Split by dataset; downsampled to the smallest
      dataset.
      PNG format
    "sd:visualPlugins":
    - image:
        tab: "Split by dataset"
        Caption: "Composition plot colored by cell cycle phase (split by dataset, downsampled)"

  umap_gr_clst_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/umap_gr_clst_spl_cnd_res_plot_png
    label: "UMAP colored by cluster (split by grouping condition, downsampled)"
    doc: |
      UMAP colored by cluster.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by group"
        Caption: "UMAP colored by cluster (split by grouping condition, downsampled)"

  cmp_gr_clst_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/cmp_gr_clst_spl_cnd_res_plot_png
    label: "Composition plot colored by cluster (split by grouping condition, downsampled)"
    doc: |
      Composition plot colored by cluster.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by group"
        Caption: "Composition plot colored by cluster (split by grouping condition, downsampled)"

  umap_gr_ph_spl_cnd_plot_png:
    type: File?
    outputSource: sc_wnn_cluster/umap_gr_ph_spl_cnd_plot_png
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
    outputSource: sc_wnn_cluster/cmp_gr_ph_spl_cnd_plot_png
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

  xpr_avg_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/xpr_avg_res_plot_png
    label: "Average gene expression"
    doc: |
      Average gene expression.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Genes of interest (expression)"
        Caption: "Average gene expression"

  xpr_dnst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/xpr_dnst_res_plot_png
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
    outputSource: sc_wnn_cluster/xpr_per_cell_plot_png
    label: "UMAP colored by gene expression (per gene)"
    doc: |
      UMAP colored by gene expression.
      All genes of interest.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Genes of interest (expression)"
        Caption: "UMAP colored by gene expression (per gene)"

  cvrg_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/cvrg_res_plot_png
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
    outputSource: sc_wnn_cluster/gse_per_cell_plot_png
    label: "UMAP colored by gene set expression score"
    doc: |
      UMAP colored by gene set expression score.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Gene sets of interest (expression)"
        Caption: "UMAP colored by gene set expression score"

  gse_avg_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/gse_avg_res_plot_png
    label: "Average gene set expression score"
    doc: |
      Average gene set expression score.
      All resolutions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Gene sets of interest (expression)"
        Caption: "Average gene set expression score"

  gse_dnst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/gse_dnst_res_plot_png
    label: "Gene set expression score density"
    doc: |
      Gene set expression score density.
      All resolutions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Gene sets of interest (expression)"
        Caption: "Gene set expression score density"

  xpr_htmp_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/xpr_htmp_res_plot_png
    label: "Gene expression heatmap (top gene markers)"
    doc: |
      Gene expression heatmap.
      Top gene markers.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Gene markers heatmap"
        Caption: "Gene expression heatmap (top gene markers)"

  xpr_htmp_res_tsv:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/xpr_htmp_res_tsv
    label: "Gene expression heatmap (top gene markers)"
    doc: |
      Gene expression heatmap.
      Top gene markers.
      TSV format.

  gene_markers_tsv:
    type: File?
    outputSource: sc_wnn_cluster/gene_markers_tsv
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
    outputSource: sc_wnn_cluster/peak_markers_tsv
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
    outputSource: sc_wnn_cluster/ucsc_cb_html_data
    label: "UCSC Cell Browser (data)"
    doc: |
      UCSC Cell Browser html data.

  ucsc_cb_html_file:
    type: File?
    outputSource: sc_wnn_cluster/ucsc_cb_html_file
    label: "UCSC Cell Browser"
    doc: |
      UCSC Cell Browser html index.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  seurat_data_rds:
    type: File
    outputSource: sc_wnn_cluster/seurat_data_rds
    label: "Seurat object in RDS format"
    doc: |
      Seurat object.
      RDS format.

  seurat_data_scope:
    type: File?
    outputSource: sc_wnn_cluster/seurat_data_scope
    label: "Seurat object in SCope compatible loom format"
    doc: |
      Seurat object.
      SCope compatible.
      Loom format.

  seurat_rna_data_cloupe:
    type: File?
    outputSource: sc_wnn_cluster/seurat_rna_data_cloupe
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
    outputSource: sc_wnn_cluster/sc_report_html_file
    label: "Analysis log"
    doc: |
      Tehcnical report.
      HTML format.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  sc_wnn_cluster_stdout_log:
    type: File
    outputSource: sc_wnn_cluster/stdout_log
    label: "Output log"
    doc: |
      Stdout log from the sc_wnn_cluster step.

  sc_wnn_cluster_stderr_log:
    type: File
    outputSource: sc_wnn_cluster/stderr_log
    label: "Error log"
    doc: |
      Stderr log from the sc_wnn_cluster step.


steps:

  sc_wnn_cluster:
    doc: |
      Clusters multiome ATAC and RNA-Seq datasets, identifies
      gene markers and differentially accessible peaks
    run: ../tools/sc-wnn-cluster.cwl
    in:
      query_data_rds: query_data_rds
      rna_dimensions: rna_dimensions
      atac_dimensions: atac_dimensions
      cluster_algorithm:
        default: "slm"
      resolution:
        source: resolution
        valueFrom: $(parse_range(self))
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
    - cell_cnts_gr_clst_res_plot_png
    - gene_umi_spl_clst_res_plot_png
    - umi_mito_spl_clst_res_plot_png
    - rna_atac_cnts_spl_clst_res_plot_png
    - tss_frgm_spl_clst_res_plot_png
    - rnadbl_gr_clst_res_plot_png
    - atacdbl_gr_clst_res_plot_png
    - vrlpdbl_gr_clst_res_plot_png
    - qc_mtrcs_dnst_gr_clst_res_plot_png
    - umap_gr_ph_spl_idnt_plot_png
    - cmp_gr_ph_spl_idnt_plot_png
    - umap_gr_ph_spl_cnd_plot_png
    - cmp_gr_ph_spl_cnd_plot_png
    - umap_gr_clst_res_plot_png
    - slh_gr_clst_res_plot_png
    - umap_gr_clst_spl_idnt_res_plot_png
    - cmp_gr_clst_spl_idnt_res_plot_png
    - umap_gr_clst_spl_ph_res_plot_png
    - cmp_gr_ph_spl_clst_res_plot_png
    - umap_gr_clst_spl_cnd_res_plot_png
    - cmp_gr_clst_spl_cnd_res_plot_png
    - gse_per_cell_plot_png
    - gse_avg_res_plot_png
    - gse_dnst_res_plot_png
    - xpr_per_cell_plot_png
    - xpr_avg_res_plot_png
    - xpr_dnst_res_plot_png
    - xpr_htmp_res_plot_png
    - cvrg_res_plot_png
    - all_plots_pdf
    - xpr_htmp_res_tsv
    - gene_markers_tsv
    - peak_markers_tsv
    - ucsc_cb_html_data
    - ucsc_cb_html_file
    - seurat_data_rds
    - seurat_rna_data_cloupe
    - seurat_data_scope
    - sc_report_html_file
    - stdout_log
    - stderr_log

  folder_pdf_plots:
    run: ../tools/files-to-folder.cwl
    in:
      input_files:
        source:
        - sc_wnn_cluster/all_plots_pdf
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

label: "Single-Cell WNN Cluster Analysis"
s:name: "Single-Cell WNN Cluster Analysis"
s:alternateName: "Single-Cell WNN Cluster Analysis"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-wnn-cluster.cwl
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
  Single-Cell WNN Cluster Analysis

  Clusters cells by similarity on the basis of both
  gene expression and chromatin accessibility data
  from the outputs of the “Single-Cell RNA-Seq
  Dimensionality Reduction Analysis” and “Single-Cell
  ATAC-Seq Dimensionality Reduction Analysis” pipelines
  run sequentially. The results of this workflow are
  used in the “Single-Cell Manual Cell Type Assignment”,
  “Single-Cell RNA-Seq Differential Expression Analysis”,
  “Single-Cell RNA-Seq Trajectory Analysis”, “Single-Cell
  Differential Abundance Analysis”,  “Single-Cell ATAC-Seq
  Differential Accessibility Analysis”, and “Single-Cell
  ATAC-Seq Genome Coverage” pipelines.