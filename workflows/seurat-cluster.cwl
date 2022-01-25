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


'sd:upstream':
  sc_rnaseq_sample:
  - "https://github.com/datirium/workflows/workflows/cellranger-aggr.cwl"
  - "https://github.com/datirium/workflows/workflows/single-cell-preprocess-cellranger.cwl"
  - "cellranger-aggr.cwl"
  - "single-cell-preprocess-cellranger.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  filtered_feature_bc_matrix_folder:
    type: File
    label: "scRNA-Seq Cellranger Experiment"
    doc: |
      Compressed folder with aggregated filtered feature-barcode matrices in MEX format
    'sd:upstreamSource': "sc_rnaseq_sample/filtered_feature_bc_matrix_folder"
    'sd:localLabel': true

  aggregation_metadata:
    type: File?
    label: "scRNA-Seq Cellranger Experiment"
    doc: |
      Aggregation metadata in CSV format.
      If not provided, we assume that upstream sc_rnaseq_sample
      was not from aggregated sample
    'sd:upstreamSource': "sc_rnaseq_sample/aggregation_metadata"
    'sd:localLabel': true

  minimum_cells:
    type: int?
    default: 5
    label: "Include genes detected in at least this many cells"
    doc: |
      Include genes detected in at least this many cells
      (applied to thoughout all datasets together).
    'sd:layout':
      advanced: true

  minimum_features:
    type: string?
    default: "250"
    label: "Include cells where at least this many genes are detected"
    doc: |
      Include cells where at least this many genes are detected.
      If multiple values provided each of them will be applied to
      the correspondent dataset.
    'sd:layout':
      advanced: true

  maximum_features:
    type: string?
    default: "5000"
    label: "Include cells with the number of genes not bigger than this value"
    doc: |
      Include cells with the number of genes not bigger than this value.
      If multiple values provided each of them will be applied to the
      correspondent dataset.
    'sd:layout':
      advanced: true

  minimum_umis:
    type: string?
    default: "500"
    label: "Include cells where at least this many UMIs are detected"
    doc: |
      Include cells where at least this many UMIs are detected.
      If multiple values provided each of them will be applied
      to the correspondent dataset.
    'sd:layout':
      advanced: true

  minimum_novelty_score:
    type: string?
    default: "0.8"
    label: "Include cells with the novelty score (the ratio of genes per cell over UMIs per cell) not lower than this value"
    doc: |
      Include cells with the novelty score (the ratio of genes per cell over UMIs per cell)
      not lower than this value (calculated as log10(genes)/log10(UMIs)). If multiple values
      provided each of them will be applied to the correspondent dataset.
    'sd:layout':
      advanced: true

  maximum_mito_perc:
    type: float?
    default: 5
    label: "Include cells with the percentage of transcripts mapped to mitochondrial genes not bigger than this value"
    doc: |
      Include cells with the percentage of transcripts mapped to mitochondrial genes not bigger than this value.
    'sd:layout':
      advanced: true

  mito_pattern:
    type: string?
    default: "^Mt-"
    label: "Pattern to identify mitochondrial genes"
    doc: |
      Pattern to identify mitochondrial genes.
    'sd:layout':
      advanced: true

  high_var_features_count:
    type: int?
    default: 3000
    label: "Number of highly variable genes to detect (used for dataset integration and dimensional reduction)"
    doc: |
      Number of highly variable genes to detect (used for dataset integration and dimensional reduction).
    'sd:layout':
      advanced: true

  dimensionality:
    type: int?
    default: 10
    label: "Number of principal components to use in UMAP projection and clustering (from 1 to 50)"
    doc: |
      Number of principal components to use in UMAP projection and clustering (from 1 to 50).
      Use Elbow plot to adjust this parameter.
    'sd:layout':
      advanced: true

  umap_spread:
    type: float?
    default: 1
    label: "Effective scale of embedded points on UMAP. Determines how clustered/clumped the embedded points are."
    doc: |
      The effective scale of embedded points on UMAP. In combination with mindist
      this determines how clustered/clumped the embedded points are.
    'sd:layout':
      advanced: true

  umap_mindist:
    type: float?
    default: 0.3
    label: "Controls how tightly the embedding is allowed compress points together on UMAP. Sensible values are in the range 0.001 to 0.5"
    doc: |
      Controls how tightly the embedding is allowed compress points together on UMAP.
      Larger values ensure embedded points are moreevenly distributed, while smaller
      values allow the algorithm to optimise more accurately with regard to local structure.
      Sensible values are in the range 0.001 to 0.5.
    'sd:layout':
      advanced: true

  umap_nneighbors:
    type: int?
    default: 30
    label: "Number of neighboring points used in UMAP. Larger values result in loss of detailed local structure."
    doc: |
      Determines the number of neighboring points used in UMAP. Larger values will result
      in more global structure being preserved at the loss of detailed local structure.
      In general this parameter should often be in the range 5 to 50.
    'sd:layout':
      advanced: true

  umap_metric:
    type:
    - "null"
    - type: enum
      symbols:
      - "euclidean"
      - "manhattan"
      - "chebyshev"
      - "minkowski"
      - "canberra"
      - "braycurtis"
      - "mahalanobis"
      - "wminkowski"
      - "seuclidean"
      - "cosine"
      - "correlation"
      - "haversine"
      - "hamming"
      - "jaccard"
      - "dice"
      - "russelrao"
      - "kulsinski"
      - "ll_dirichlet"
      - "hellinger"
      - "rogerstanimoto"
      - "sokalmichener"
      - "sokalsneath"
      - "yule"
    default: "cosine"
    label: "The metric to use to compute distances in high dimensional space for UMAP"
    doc: |
      The metric to use to compute distances in high dimensional space for UMAP.
    'sd:layout':
      advanced: true

  umap_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "uwot"
      - "uwot-learn"
      - "umap-learn"
    default: "uwot"
    label: "UMAP implementation to run"
    doc: |
      UMAP implementation to run.
    'sd:layout':
      advanced: true

  cluster_metric:
    type:
    - "null"
    - type: enum
      symbols:
      - "euclidean"
      - "cosine"
      - "manhattan"
      - "hamming"
    default: "euclidean"
    label: "Distance metric used by the nearest neighbors algorithm when running clustering"
    doc: |
      Distance metric used by the nearest neighbors algorithm when running clustering.
    'sd:layout':
      advanced: true

  resolution:
    type: string?
    default: "0.1"
    label: "Comma or space separated list of clustering resolutions"
    doc: |
      Comma or space separated list of clustering resolutions
    'sd:layout':
      advanced: true

  minimum_logfc:
    type: float?
    default: 0.25
    label: "Include only those genes that on average have log fold change difference in expression between every tested pair of clusters not lower than this value"
    doc: |
      Include only those genes that on average have log fold change difference in
      expression between every tested pair of clusters not lower than this value.
    'sd:layout':
      advanced: true

  minimum_pct:
    type: float?
    default: 0.1
    label: "Include only those genes that are detected in not lower than this fraction of cells in either of the two tested clusters"
    doc: |
      Include only those genes that are detected in not lower than
      this fraction of cells in either of the two tested clusters.
    'sd:layout':
      advanced: true

  test_use:
    type:
    - "null"
    - type: enum
      symbols:
      - "wilcox"
      - "bimod"
      - "roc"
      - "t"
      - "negbinom"
      - "poisson"
      - "LR"
      - "MAST"
      - "DESeq2"
    default: "wilcox"
    label: "Statistical test to use for gene markers identification"
    doc: |
      Statistical test to use for gene markers identification.
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 6
    label: "Threads number to use"
    doc: |
      Threads number
    'sd:layout':
      advanced: true

  species:
    type:
    - "null"
    - type: enum
      symbols:
      - "hs"
      - "mm"
      - "none"
    default: "none"
    label: "Species for gene name conversion when running cell type prediction"
    doc: |
      Select species for gene name conversion when running cell type prediction
      with Garnett classifier.
      If "none" - do not convert gene names
    'sd:layout':
      advanced: true

  regress_cellcycle:
    type: boolean?
    default: false
    label: "Regress cell cycle as a confounding source of variation"
    doc: |
      Regress cell cycle as a confounding source of variation.
    'sd:layout':
      advanced: true
      
  regress_mito_perc:
    type: boolean?
    default: false
    label: "Regress mitochondrial gene expression as a confounding source of variation"
    doc: |
      Regress mitochondrial gene expression as a confounding source
      of variation.
    'sd:layout':
      advanced: true

  only_positive_markers:
    type: boolean?
    default: false
    label: "Report only positive gene markers"
    doc: |
      Report only positive gene markers.
    'sd:layout':
      advanced: true

  no_sct:
    type: boolean?
    default: false
    label: "Use LogNormalize instead of SCTransform when integrating datasets"
    doc: |
      Do not use SCTransform when running datasets integration. Use LogNormalize instead.
    'sd:layout':
      advanced: true

  selected_features:
    type: string?
    default: null
    label: "Comma or space separated list of genes of interest"
    doc: |
      Comma or space separated list of genes of interest.
      Default: do not highlight any features
    'sd:layout':
      advanced: true

  conditions_data:
    type: File?
    label: "TSV/CSV file to define datasets conditions with 'library_id' and 'condition' columns. Rows order should correspond to the aggregation metadata."
    doc: |
      Path to the TSV/CSV file to define datasets grouping. First column -
      'library_id' with the values provided in the same order as in the
      correspondent column of the --identity file, second column 'condition'.
      If not provided, each dataset is assigned to its own
      biological condition

  barcodes_data:
    type: File?
    label: "Headerless TSV/CSV file with cell barcodes (one barcode per line) to prefilter input data"
    doc: |
      Path to the headerless TSV/CSV file with selected barcodes
      (one per line) to prefilter input feature-barcode matrices.
      If not provided, use all cells
    'sd:layout':
      advanced: true

  cell_cycle_data:
    type: File?
    label: "TSV/CSV file with cell cycle data with 'phase' and 'gene_id' columns"
    doc: |
      TSV/CSV file with cell cycle data. First column - 'phase', second column 'gene_id'.
      If not provided, skip cell cycle score assignment
    'sd:layout':
      advanced: true

  classifier_rds:
    type: File?
    label: "Garnett classifier rds file for cell type prediction"
    doc: |
      Path to the Garnett classifier rds file for cell type prediction.
      If not provided, skip cell type prediction
    'sd:layout':
      advanced: true


outputs:

  raw_cell_count_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_cell_count_plot_png
    label: "Number of cells per dataset (not filtered)"
    doc: |
      Number of cells per dataset (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Number of cells per dataset (not filtered)'

  raw_cell_count_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_cell_count_plot_pdf
    label: "Number of cells per dataset (not filtered)"
    doc: |
      Number of cells per dataset (not filtered).
      PDF format


  raw_umi_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_umi_dnst_spl_by_cond_plot_png
    label: "Split by condition UMI density per cell (not filtered)"
    doc: |
      Split by condition UMI density per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Split by condition UMI density per cell (not filtered)'

  raw_umi_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_umi_dnst_spl_by_cond_plot_pdf
    label: "Split by condition UMI density per cell (not filtered)"
    doc: |
      Split by condition UMI density per cell (not filtered).
      PDF format


  raw_gene_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_gene_dnst_spl_by_cond_plot_png
    label: "Split by condition gene density per cell (not filtered)"
    doc: |
      Split by condition gene density per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Split by condition gene density per cell (not filtered)'

  raw_gene_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_gene_dnst_spl_by_cond_plot_pdf
    label: "Split by condition gene density per cell (not filtered)"
    doc: |
      Split by condition gene density per cell (not filtered).
      PDF format


  raw_gene_umi_corr_spl_by_ident_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_gene_umi_corr_spl_by_ident_plot_png
    label: "Split by identity genes vs UMIs per cell correlation (not filtered)"
    doc: |
      Split by identity genes vs UMIs per cell correlation (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Split by identity genes vs UMIs per cell correlation (not filtered)'

  raw_gene_umi_corr_spl_by_ident_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_gene_umi_corr_spl_by_ident_plot_pdf
    label: "Split by identity genes vs UMIs per cell correlation (not filtered)"
    doc: |
      Split by identity genes vs UMIs per cell correlation (not filtered).
      PDF format


  raw_mito_perc_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_mito_perc_dnst_spl_by_cond_plot_png
    label: "Split by condition density of transcripts mapped to mitochondrial genes per cell (not filtered)"
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Split by condition density of transcripts mapped to mitochondrial genes per cell (not filtered)'

  raw_mito_perc_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_mito_perc_dnst_spl_by_cond_plot_pdf
    label: "Split by condition density of transcripts mapped to mitochondrial genes per cell (not filtered)"
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (not filtered).
      PDF format


  raw_nvlt_score_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_nvlt_score_dnst_spl_by_cond_plot_png
    label: "Split by condition novelty score density per cell (not filtered)"
    doc: |
      Split by condition novelty score density per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Split by condition novelty score density per cell (not filtered)'

  raw_nvlt_score_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_nvlt_score_dnst_spl_by_cond_plot_pdf
    label: "Split by condition novelty score density per cell (not filtered)"
    doc: |
      Split by condition novelty score density per cell (not filtered).
      PDF format


  raw_qc_mtrcs_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_qc_mtrcs_plot_png
    label: "QC metrics densities per cell (not filtered)"
    doc: |
      QC metrics densities per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'QC metrics densities per cell (not filtered)'

  raw_qc_mtrcs_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_qc_mtrcs_plot_pdf
    label: "QC metrics densities per cell (not filtered)"
    doc: |
      QC metrics densities per cell (not filtered).
      PDF format

  raw_miqc_mtrcs_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_miqc_mtrcs_plot_png
    label: "MiQC prediction of the compromised cells level (not filtered)"
    doc: |
      MiQC prediction of the compromised cells level (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'MiQC prediction of the compromised cells level (not filtered)'

  raw_miqc_mtrcs_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_miqc_mtrcs_plot_pdf
    label: "MiQC prediction of the compromised cells level (not filtered)"
    doc: |
      MiQC prediction of the compromised cells level (not filtered).
      PDF format

  raw_qc_mtrcs_gr_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_qc_mtrcs_gr_by_cond_plot_png
    label: "Grouped by condition QC metrics densities per cell (not filtered)"
    doc: |
      Grouped by condition QC metrics densities per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Grouped by condition QC metrics densities per cell (not filtered)'

  raw_qc_mtrcs_gr_by_cond_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_qc_mtrcs_gr_by_cond_plot_pdf
    label: "Grouped by condition QC metrics densities per cell (not filtered)"
    doc: |
      Grouped by condition QC metrics densities per cell (not filtered).
      PDF format


  fltr_cell_count_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_cell_count_plot_png
    label: "Number of cells per dataset (filtered)"
    doc: |
      Number of cells per dataset (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Number of cells per dataset (filtered)'

  fltr_cell_count_plot_pdf:
    type: File?
    outputSource: seurat_cluster/fltr_cell_count_plot_pdf
    label: "Number of cells per dataset (filtered)"
    doc: |
      Number of cells per dataset (filtered).
      PDF format


  fltr_umi_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_umi_dnst_spl_by_cond_plot_png
    label: "Split by condition UMI density per cell (filtered)"
    doc: |
      Split by condition UMI density per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Split by condition UMI density per cell (filtered)'

  fltr_umi_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputSource: seurat_cluster/fltr_umi_dnst_spl_by_cond_plot_pdf
    label: "Split by condition UMI density per cell (filtered)"
    doc: |
      Split by condition UMI density per cell (filtered).
      PDF format


  fltr_gene_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_gene_dnst_spl_by_cond_plot_png
    label: "Split by condition gene density per cell (filtered)"
    doc: |
      Split by condition gene density per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Split by condition gene density per cell (filtered)'

  fltr_gene_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputSource: seurat_cluster/fltr_gene_dnst_spl_by_cond_plot_pdf
    label: "Split by condition gene density per cell (filtered)"
    doc: |
      Split by condition gene density per cell (filtered).
      PDF format


  fltr_gene_umi_corr_spl_by_ident_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_gene_umi_corr_spl_by_ident_plot_png
    label: "Split by identity genes vs UMIs per cell correlation (filtered)"
    doc: |
      Split by identity genes vs UMIs per cell correlation (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Split by identity genes vs UMIs per cell correlation (filtered)'

  fltr_gene_umi_corr_spl_by_ident_plot_pdf:
    type: File?
    outputSource: seurat_cluster/fltr_gene_umi_corr_spl_by_ident_plot_pdf
    label: "Split by identity genes vs UMIs per cell correlation (filtered)"
    doc: |
      Split by identity genes vs UMIs per cell correlation (filtered).
      PDF format


  fltr_mito_perc_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_mito_perc_dnst_spl_by_cond_plot_png
    label: "Split by condition density of transcripts mapped to mitochondrial genes per cell (filtered)"
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Split by condition density of transcripts mapped to mitochondrial genes per cell (filtered)'

  fltr_mito_perc_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputSource: seurat_cluster/fltr_mito_perc_dnst_spl_by_cond_plot_pdf
    label: "Split by condition density of transcripts mapped to mitochondrial genes per cell (filtered)"
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (filtered).
      PDF format


  fltr_nvlt_score_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_nvlt_score_dnst_spl_by_cond_plot_png
    label: "Split by condition novelty score density per cell (filtered)"
    doc: |
      Split by condition novelty score density per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Split by condition novelty score density per cell (filtered)'

  fltr_nvlt_score_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputSource: seurat_cluster/fltr_nvlt_score_dnst_spl_by_cond_plot_pdf
    label: "Split by condition novelty score density per cell (filtered)"
    doc: |
      Split by condition novelty score density per cell (filtered).
      PDF format


  fltr_qc_mtrcs_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_qc_mtrcs_plot_png
    label: "QC metrics densities per cell (filtered)"
    doc: |
      QC metrics densities per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'QC metrics densities per cell (filtered)'

  fltr_qc_mtrcs_plot_pdf:
    type: File?
    outputSource: seurat_cluster/fltr_qc_mtrcs_plot_pdf
    label: "QC metrics densities per cell (filtered)"
    doc: |
      QC metrics densities per cell (filtered).
      PDF format


  fltr_qc_mtrcs_gr_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_qc_mtrcs_gr_by_cond_plot_png
    label: "Grouped by condition QC metrics densities per cell (filtered)"
    doc: |
      Grouped by condition QC metrics densities per cell (filtered).
      PDF format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Grouped by condition QC metrics densities per cell (filtered)'

  fltr_qc_mtrcs_gr_by_cond_plot_pdf:
    type: File?
    outputSource: seurat_cluster/fltr_qc_mtrcs_gr_by_cond_plot_pdf
    label: "Grouped by condition QC metrics densities per cell (filtered)"
    doc: |
      Grouped by condition QC metrics densities per cell (filtered).
      PDF format


  fltr_pca_spl_by_ph_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_pca_spl_by_ph_plot_png
    label: "Split by cell cycle phase PCA of filtered unintegrated/scaled datasets"
    doc: |
      Split by cell cycle phase PCA of filtered unintegrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Split by cell cycle phase PCA of filtered unintegrated/scaled datasets'

  fltr_pca_spl_by_ph_plot_pdf:
    type: File?
    outputSource: seurat_cluster/fltr_pca_spl_by_ph_plot_pdf
    label: "Split by cell cycle phase PCA of filtered unintegrated/scaled datasets"
    doc: |
      Split by cell cycle phase PCA of filtered unintegrated/scaled datasets.
      PDF format


  fltr_pca_spl_by_mito_perc_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_pca_spl_by_mito_perc_plot_png
    label: "Split by level of transcripts mapped to mitochondrial genes PCA of filtered unintegrated/scaled datasets"
    doc: |
      Split by level of transcripts mapped to mitochondrial genes PCA of filtered unintegrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Split by level of transcripts mapped to mitochondrial genes PCA of filtered unintegrated/scaled datasets'

  fltr_pca_spl_by_mito_perc_plot_pdf:
    type: File?
    outputSource: seurat_cluster/fltr_pca_spl_by_mito_perc_plot_pdf
    label: "Split by level of transcripts mapped to mitochondrial genes PCA of filtered unintegrated/scaled datasets"
    doc: |
      Split by level of transcripts mapped to mitochondrial genes PCA of filtered unintegrated/scaled datasets.
      PDF format


  fltr_umap_spl_by_idnt_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_umap_spl_by_idnt_plot_png
    label: "Split by identity UMAP projected PCA of filtered unintegrated/scaled datasets"
    doc: |
      Split by identity UMAP projected PCA of filtered unintegrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Split by identity UMAP projected PCA of filtered unintegrated/scaled datasets'

  fltr_umap_spl_by_idnt_plot_pdf:
    type: File?
    outputSource: seurat_cluster/fltr_umap_spl_by_idnt_plot_pdf
    label: "Split by identity UMAP projected PCA of filtered unintegrated/scaled datasets"
    doc: |
      Split by identity UMAP projected PCA of filtered unintegrated/scaled datasets.
      PDF format


  ntgr_elbow_plot_png:
    type: File?
    outputSource: seurat_cluster/ntgr_elbow_plot_png
    label: "Elbow plot from PCA of filtered integrated/scaled datasets"
    doc: |
      Elbow plot from PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 3. Dimensionality evaluation'
        Caption: 'Elbow plot from PCA of filtered integrated/scaled datasets'

  ntgr_elbow_plot_pdf:
    type: File?
    outputSource: seurat_cluster/ntgr_elbow_plot_pdf
    label: "Elbow plot from PCA of filtered integrated/scaled datasets"
    doc: |
      Elbow plot from PCA of filtered integrated/scaled datasets.
      PDF format


  ntgr_pca_plot_png:
    type: File?
    outputSource: seurat_cluster/ntgr_pca_plot_png
    label: "PCA of filtered integrated/scaled datasets"
    doc: |
      PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 3. Dimensionality evaluation'
        Caption: 'PCA of filtered integrated/scaled datasets'

  ntgr_pca_plot_pdf:
    type: File?
    outputSource: seurat_cluster/ntgr_pca_plot_pdf
    label: "PCA of filtered integrated/scaled datasets"
    doc: |
      PCA of filtered integrated/scaled datasets.
      PDF format


  ntgr_pca_heatmap_png:
    type: File?
    outputSource: seurat_cluster/ntgr_pca_heatmap_png
    label: "Genes per cells expression heatmap sorted by their PC scores from PCA of filtered integrated/scaled datasets"
    doc: |
      Genes per cells expression heatmap sorted by their PC scores from PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 3. Dimensionality evaluation'
        Caption: 'Genes per cells expression heatmap sorted by their PC scores from PCA of filtered integrated/scaled datasets'

  ntgr_pca_heatmap_pdf:
    type: File?
    outputSource: seurat_cluster/ntgr_pca_heatmap_pdf
    label: "Genes per cells expression heatmap sorted by their PC scores from PCA of filtered integrated/scaled datasets"
    doc: |
      Genes per cells expression heatmap sorted by their PC scores from PCA of filtered integrated/scaled datasets.
      PDF format


  ntgr_pca_loadings_plot_png:
    type: File?
    outputSource: seurat_cluster/ntgr_pca_loadings_plot_png
    label: "PC scores of the most variant genes from PCA of filtered integrated/scaled datasets"
    doc: |
      PC scores of the most variant genes from PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 3. Dimensionality evaluation'
        Caption: 'PC scores of the most variant genes from PCA of filtered integrated/scaled datasets'

  ntgr_pca_loadings_plot_pdf:
    type: File?
    outputSource: seurat_cluster/ntgr_pca_loadings_plot_pdf
    label: "PC scores of the most variant genes from PCA of filtered integrated/scaled datasets"
    doc: |
      PC scores of the most variant genes from PCA of filtered integrated/scaled datasets.
      PDF format


  ntgr_umap_spl_by_idnt_plot_png:
    type: File?
    outputSource: seurat_cluster/ntgr_umap_spl_by_idnt_plot_png
    label: "Split by identity UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      Split by identity UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 4. Integration QC'
        Caption: 'Split by identity UMAP projected PCA of filtered integrated/scaled datasets'

  ntgr_umap_spl_by_idnt_plot_pdf:
    type: File?
    outputSource: seurat_cluster/ntgr_umap_spl_by_idnt_plot_pdf
    label: "Split by identity UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      Split by identity UMAP projected PCA of filtered integrated/scaled datasets.
      PDF format


  clst_umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/clst_umap_res_plot_png
    label: "Clustered UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      Clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 5. Clustering'
        Caption: 'Clustered UMAP projected PCA of filtered integrated/scaled datasets'

  clst_umap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/clst_umap_res_plot_pdf
    label: "Clustered UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      Clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PDF format


  clst_umap_spl_by_cond_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/clst_umap_spl_by_cond_res_plot_png
    label: "Split by condition clustered UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      Split by condition clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 5. Clustering'
        Caption: 'Split by condition clustered UMAP projected PCA of filtered integrated/scaled datasets'

  clst_umap_spl_by_cond_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/clst_umap_spl_by_cond_res_plot_pdf
    label: "Split by condition clustered UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      Split by condition clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PDF format


  clst_umap_ctype_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/clst_umap_ctype_res_plot_png
    label: "Grouped by predicted cell types UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      Grouped by predicted cell types UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 5. Clustering'
        Caption: 'Grouped by predicted cell types UMAP projected PCA of filtered integrated/scaled datasets'

  clst_umap_ctype_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/clst_umap_ctype_res_plot_pdf
    label: "Grouped by predicted cell types UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      Grouped by predicted cell types UMAP projected PCA of filtered integrated/scaled datasets.
      PDF format


  clst_umap_spl_by_ph_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/clst_umap_spl_by_ph_res_plot_png
    label: "Split by cell cycle phase clustered UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      Split by cell cycle phase clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 4. Integration QC'
        Caption: 'Split by cell cycle phase clustered UMAP projected PCA of filtered integrated/scaled datasets'

  clst_umap_spl_by_ph_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/clst_umap_spl_by_ph_res_plot_pdf
    label: "Split by cell cycle phase clustered UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      Split by cell cycle phase clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PDF format


  clst_qc_mtrcs_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/clst_qc_mtrcs_res_plot_png
    label: "QC metrics for clustered UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      QC metrics for clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 4. Integration QC'
        Caption: 'QC metrics for clustered UMAP projected PCA of filtered integrated/scaled datasets'

  clst_qc_mtrcs_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/clst_qc_mtrcs_res_plot_pdf
    label: "QC metrics for clustered UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      QC metrics for clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PDF format


  expr_avg_per_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_avg_per_clst_res_plot_png
    label: "Scaled average log normalized gene expression per cluster of filtered integrated/scaled datasets"
    doc: |
      Scaled average log normalized gene expression per cluster of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 6. Gene expression'
        Caption: 'Scaled average log normalized gene expression per cluster of filtered integrated/scaled datasets'

  expr_avg_per_clst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_avg_per_clst_res_plot_pdf
    label: "Scaled average log normalized gene expression per cluster of filtered integrated/scaled datasets"
    doc: |
      Scaled average log normalized gene expression per cluster of filtered integrated/scaled datasets.
      PDF format


  expr_per_clst_cell_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_per_clst_cell_res_plot_png
    label: "Log normalized gene expression per cell of clustered filtered integrated/scaled datasets"
    doc: |
      Log normalized gene expression per cell of clustered filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 6. Gene expression'
        Caption: 'Log normalized gene expression per cell of clustered filtered integrated/scaled datasets'

  expr_per_clst_cell_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_per_clst_cell_res_plot_pdf
    label: "Log normalized gene expression per cell of clustered filtered integrated/scaled datasets"
    doc: |
      Log normalized gene expression per cell of clustered filtered integrated/scaled datasets.
      PDF format


  expr_clst_heatmap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_clst_heatmap_res_plot_png
    label: "Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets"
    doc: |
      Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 6. Gene expression'
        Caption: 'Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets'

  expr_clst_heatmap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_clst_heatmap_res_plot_pdf
    label: "Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets"
    doc: |
      Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets.
      PDF format


  expr_dnst_per_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_dnst_per_clst_res_plot_png
    label: "Log normalized gene expression densities per cluster of filtered integrated/scaled datasets"
    doc: |
      Log normalized gene expression densities per cluster of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 6. Gene expression'
        Caption: 'Log normalized gene expression densities per cluster of filtered integrated/scaled datasets'

  expr_dnst_per_clst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_dnst_per_clst_res_plot_pdf
    label: "Log normalized gene expression densities per cluster of filtered integrated/scaled datasets"
    doc: |
      Log normalized gene expression densities per cluster of filtered integrated/scaled datasets.
      PDF format


  expr_avg_per_ctype_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_avg_per_ctype_res_plot_png
    label: "Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets"
    doc: |
      Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 6. Gene expression'
        Caption: 'Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets'

  expr_avg_per_ctype_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_avg_per_ctype_res_plot_pdf
    label: "Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets"
    doc: |
      Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets.
      PDF format


  expr_per_ctype_cell_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_per_ctype_cell_res_plot_png
    label: "Log normalized gene expression per cell of clustered filtered integrated/scaled datasets with predicted cell types"
    doc: |
      Log normalized gene expression per cell of clustered filtered/scaled integrated datasets with predicted cell types.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 6. Gene expression'
        Caption: 'Log normalized gene expression per cell of clustered filtered/scaled integrated datasets with predicted cell types'

  expr_per_ctype_cell_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_per_ctype_cell_res_plot_pdf
    label: "Log normalized gene expression per cell of clustered filtered integrated/scaled datasets with predicted cell types"
    doc: |
      Log normalized gene expression per cell of clustered filtered integrated/scaled datasets with predicted cell types.
      PDF format


  expr_ctype_heatmap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_ctype_heatmap_res_plot_png
    label: "Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets with predicted cell types"
    doc: |
      Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets with predicted cell types.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 6. Gene expression'
        Caption: 'Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets with predicted cell types'

  expr_ctype_heatmap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_ctype_heatmap_res_plot_pdf
    label: "Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets with predicted cell types"
    doc: |
      Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets with predicted cell types.
      PDF format


  expr_dnst_per_ctype_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_dnst_per_ctype_res_plot_png
    label: "Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets"
    doc: |
      Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 6. Gene expression'
        Caption: 'Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets'

  expr_dnst_per_ctype_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_dnst_per_ctype_res_plot_pdf
    label: "Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets"
    doc: |
      Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets.
      PDF format


  seurat_clst_data_rds:
    type: File
    outputSource: seurat_cluster/seurat_clst_data_rds
    label: "Clustered filtered integrated/scaled Seurat data"
    doc: |
      Clustered filtered integrated Seurat data.
      RDS format


  clst_pttv_gene_markers:
    type: File
    outputSource: seurat_cluster/clst_pttv_gene_markers
    label: "Putative gene markers file for all clusters and all resolutions"
    doc: |
      Putative gene markers file for all clusters and all resolutions.
      TSV format
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Step 7. Putative gene markers'
        Title: 'Putative gene markers'


  clst_csrvd_gene_markers:
    type: File
    outputSource: seurat_cluster/clst_csrvd_gene_markers
    label: "Conserved gene markers file for all clusters and all resolutions"
    doc: |
      Conserved gene markers file for all clusters and all resolutions.
      TSV format
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Step 8. Conserved gene markers'
        Title: 'Conserved gene markers'


  compressed_cellbrowser_config_data:
    type: File
    outputSource: compress_cellbrowser_config_data/compressed_folder
    label: "Compressed directory with UCSC Cellbrowser configuration data"
    doc: |
      Compressed directory with UCSC Cellbrowser configuration data

  cellbrowser_html_data:
    type: Directory
    outputSource: seurat_cluster/cellbrowser_html_data
    label: "Directory with UCSC Cellbrowser formatted html data"
    doc: |
      Directory with UCSC Cellbrowser formatted html data

  cellbrowser_html_file:
    type: File
    outputSource: seurat_cluster/cellbrowser_html_file
    label: "Open in UCSC Cell Browser"
    doc: |
      HTML index file from the directory with UCSC Cellbrowser formatted html data
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"


  seurat_cluster_stdout_log:
    type: File
    outputSource: seurat_cluster/stdout_log
    label: stdout log generated by Seurat
    doc: |
      stdout log generated by Seurat

  seurat_cluster_stderr_log:
    type: File
    outputSource: seurat_cluster/stderr_log
    label: stderr log generated by Seurat
    doc: |
      stderr log generated by Seurat


steps:

  uncompress_feature_bc_matrices:
    in:
      compressed: filtered_feature_bc_matrix_folder
    out:
    - uncompressed
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap:v0.0.3
      inputs:
        compressed:
          type: File
          inputBinding:
            position: 1
      outputs:
        uncompressed:
          type: Directory
          outputBinding:
            glob: "*"
      baseCommand: ["tar", "xzf"]

  seurat_cluster:
    run: ../tools/seurat-cluster.cwl
    in:
      feature_bc_matrices_folder: uncompress_feature_bc_matrices/uncompressed
      aggregation_metadata: aggregation_metadata
      cell_cycle_data: cell_cycle_data
      conditions_data: conditions_data
      classifier_rds: classifier_rds
      species: species
      barcodes_data: barcodes_data
      minimum_cells: minimum_cells
      minimum_features:
        source: minimum_features
        valueFrom: $(split_numbers(self))
      maximum_features:
        source: maximum_features
        valueFrom: $(split_numbers(self))
      selected_features:
        source: selected_features
        valueFrom: $(split_features(self))
      minimum_umis:
        source: minimum_umis
        valueFrom: $(split_numbers(self))
      minimum_novelty_score:
        source: minimum_novelty_score
        valueFrom: $(split_numbers(self))
      maximum_mito_perc: maximum_mito_perc
      mito_pattern: mito_pattern
      regress_cellcycle: regress_cellcycle
      regress_mito_perc: regress_mito_perc
      high_var_features_count: high_var_features_count
      dimensionality: dimensionality
      umap_spread: umap_spread
      umap_mindist: umap_mindist
      umap_nneighbors: umap_nneighbors
      umap_metric: umap_metric
      umap_method: umap_method
      no_sct: no_sct
      cluster_metric: cluster_metric
      resolution:
        source: resolution
        valueFrom: $(split_numbers(self))
      minimum_logfc: minimum_logfc
      minimum_pct: minimum_pct
      only_positive_markers: only_positive_markers
      test_use: test_use
      export_pdf_plots:
        default: true
      export_rds_data:
        default: true
      threads: threads
    out:
    - raw_cell_count_plot_png
    - raw_cell_count_plot_pdf
    - raw_umi_dnst_spl_by_cond_plot_png
    - raw_umi_dnst_spl_by_cond_plot_pdf
    - raw_gene_dnst_spl_by_cond_plot_png
    - raw_gene_dnst_spl_by_cond_plot_pdf
    - raw_gene_umi_corr_spl_by_ident_plot_png
    - raw_gene_umi_corr_spl_by_ident_plot_pdf
    - raw_mito_perc_dnst_spl_by_cond_plot_png
    - raw_mito_perc_dnst_spl_by_cond_plot_pdf
    - raw_nvlt_score_dnst_spl_by_cond_plot_png
    - raw_nvlt_score_dnst_spl_by_cond_plot_pdf
    - raw_qc_mtrcs_plot_png
    - raw_qc_mtrcs_plot_pdf
    - raw_miqc_mtrcs_plot_png
    - raw_miqc_mtrcs_plot_pdf
    - raw_qc_mtrcs_gr_by_cond_plot_png
    - raw_qc_mtrcs_gr_by_cond_plot_pdf
    - fltr_cell_count_plot_png
    - fltr_cell_count_plot_pdf
    - fltr_umi_dnst_spl_by_cond_plot_png
    - fltr_umi_dnst_spl_by_cond_plot_pdf
    - fltr_gene_dnst_spl_by_cond_plot_png
    - fltr_gene_dnst_spl_by_cond_plot_pdf
    - fltr_gene_umi_corr_spl_by_ident_plot_png
    - fltr_gene_umi_corr_spl_by_ident_plot_pdf
    - fltr_mito_perc_dnst_spl_by_cond_plot_png
    - fltr_mito_perc_dnst_spl_by_cond_plot_pdf
    - fltr_nvlt_score_dnst_spl_by_cond_plot_png
    - fltr_nvlt_score_dnst_spl_by_cond_plot_pdf
    - fltr_qc_mtrcs_plot_png
    - fltr_qc_mtrcs_plot_pdf
    - fltr_qc_mtrcs_gr_by_cond_plot_png
    - fltr_qc_mtrcs_gr_by_cond_plot_pdf
    - fltr_pca_spl_by_ph_plot_png
    - fltr_pca_spl_by_ph_plot_pdf
    - fltr_pca_spl_by_mito_perc_plot_png
    - fltr_pca_spl_by_mito_perc_plot_pdf
    - fltr_umap_spl_by_idnt_plot_png
    - fltr_umap_spl_by_idnt_plot_pdf
    - ntgr_elbow_plot_png
    - ntgr_elbow_plot_pdf
    - ntgr_pca_plot_png
    - ntgr_pca_plot_pdf
    - ntgr_pca_heatmap_png
    - ntgr_pca_heatmap_pdf
    - ntgr_pca_loadings_plot_png
    - ntgr_pca_loadings_plot_pdf
    - ntgr_umap_spl_by_idnt_plot_png
    - ntgr_umap_spl_by_idnt_plot_pdf
    - clst_umap_res_plot_png
    - clst_umap_res_plot_pdf
    - clst_umap_spl_by_cond_res_plot_png
    - clst_umap_spl_by_cond_res_plot_pdf
    - clst_umap_ctype_res_plot_png
    - clst_umap_ctype_res_plot_pdf
    - clst_umap_spl_by_ph_res_plot_png
    - clst_umap_spl_by_ph_res_plot_pdf
    - clst_qc_mtrcs_res_plot_png
    - clst_qc_mtrcs_res_plot_pdf
    - clst_pttv_gene_markers   
    - clst_csrvd_gene_markers
    - expr_avg_per_clst_res_plot_png
    - expr_avg_per_clst_res_plot_pdf
    - expr_per_clst_cell_res_plot_png
    - expr_per_clst_cell_res_plot_pdf
    - expr_clst_heatmap_res_plot_png
    - expr_clst_heatmap_res_plot_pdf
    - expr_dnst_per_clst_res_plot_png
    - expr_dnst_per_clst_res_plot_pdf
    - expr_avg_per_ctype_res_plot_png
    - expr_avg_per_ctype_res_plot_pdf
    - expr_per_ctype_cell_res_plot_png
    - expr_per_ctype_cell_res_plot_pdf
    - expr_ctype_heatmap_res_plot_png
    - expr_ctype_heatmap_res_plot_pdf
    - expr_dnst_per_ctype_res_plot_png
    - expr_dnst_per_ctype_res_plot_pdf
    - seurat_clst_data_rds
    - cellbrowser_config_data
    - cellbrowser_html_data
    - cellbrowser_html_file
    - stdout_log
    - stderr_log

  compress_cellbrowser_config_data:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: seurat_cluster/cellbrowser_config_data
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Seurat Cluster"
label: "Seurat Cluster"
s:alternateName: "Runs filtering, integration, and clustering analyses for Cell Ranger Count Gene Expression or Cell Ranger Aggregate experiments"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/seurat-cluster.cwl
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
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  Seurat Cluster
  ==============

  Runs filtering, integration, and clustering analyses for Cell Ranger
  Count Gene Expression or Cell Ranger Aggregate experiments.