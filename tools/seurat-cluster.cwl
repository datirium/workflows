cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/seurat:v0.0.12


inputs:

  feature_bc_matrices_folder:
    type:
    - Directory
    - type: array
      items: Directory
    inputBinding:
      prefix: "--mex"
    doc: |
      Path to the folder with not normalized aggregated feature-barcode matrix
      from Cell Ranger Aggregate in MEX format. If multiple locations provided
      data is assumed to be not aggregated (outputs from multiple Cell Ranger
      Count runs) and will be merged.

  aggregation_metadata:
    type: File
    inputBinding:
      prefix: "--identity"
    doc: |
      Path to the metadata TSV/CSV file to set the datasets identities.
      If --mex points to the Cell Ranger Aggregate outputs, the aggregation.csv
      file can be used as well. If multiple locations were provided through --mex,
      the file should include at least one column - 'library_id', and be sorted
      based on the the order of locations provided in --mex.

  conditions_data:
    type: File?
    inputBinding:
      prefix: "--condition"
    doc: |
      Path to the TSV/CSV file to define datasets grouping. First column -
      'library_id' with the values provided in the same order as in the
      correspondent column of the --identity file, second column 'condition'.
      Default: each dataset is assigned to a separate group.

  classifier_rds:
    type: File?
    inputBinding:
      prefix: "--classifier"
    doc: |
      Path to the Garnett classifier RDS file for cell type prediction.
      Default: skip cell type prediction.

  cell_cycle_data:
    type: File?
    inputBinding:
      prefix: "--cellcycle"
    doc: |
      Path to the TSV/CSV file with cell cycle data. First column - 'phase',
      second column 'gene_id'. Default: skip cell cycle score assignment.

  barcodes_data:
    type: File?
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the headerless TSV/CSV file with the list of barcodes to select
      cells of interest (one barcode per line). Prefilters input feature-barcode
      matrix to include only selected cells. Default: use all cells.

  minimum_cells:
    type: int?
    inputBinding:
      prefix: "--mincells"
    doc: |
      Include only features detected in at least this many cells. Applied to
      aggregated feature-barcode matrix from Cell Ranger Aggregate. Ignored
      when --mex points to the locations of multiple Cell Ranger Count runs.
      Default: 5

  minimum_features:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--minfeatures"
    doc: |
      Include cells where at least this many features are detected. If multiple
      values provided each of them will be applied to the correspondent dataset
      from the --mex input.
      Default: 250 (applied to all datasets)

  maximum_features:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--maxfeatures"
    doc: |
      Include cells with the number of features not bigger than this value. If
      multiple values provided each of them will be applied to the correspondent
      dataset from the --mex input.
      Default: 5000 (applied to all datasets)

  minimum_umis:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--minumi"
    doc: |
      Include cells where at least this many UMIs (transcripts) are detected. If
      multiple values provided each of them will be applied to the correspondent
      dataset from the --mex input.
      Default: 500 (applied to all datasets)

  minimum_novelty_score:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--minnovelty"
    doc: |
      Include cells with the novelty score not lower than this value, calculated as
      log10(genes)/log10(UMIs). If multiple values provided each of them will be
      applied to the correspondent dataset from the --mex input.
      Default: 0.8 (applied to all datasets)

  maximum_mito_perc:
    type: float?
    inputBinding:
      prefix: "--maxmt"
    doc: |
      Include cells with the percentage of transcripts mapped to mitochondrial genes
      not bigger than this value.
      Default: 5

  mito_pattern:
    type: string?
    inputBinding:
      prefix: "--mitopattern"
    doc: |
      Regex pattern to identify mitochondrial genes.
      Default: '^Mt-'

  selected_features:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--features"
    doc: |
      Features of interest to evaluate expression.
      Default: None

  regress_cellcycle:
    type: boolean?
    inputBinding:
      prefix: "--regresscellcycle"
    doc: |
      Regress cell cycle as a confounding source of variation.
      Default: false

  regress_mito_perc:
    type: boolean?
    inputBinding:
      prefix: "--regressmt"
    doc: |
      Regress mitochondrial genes expression as a confounding source of variation.
      Default: false

  high_var_features_count:
    type: int?
    inputBinding:
      prefix: "--highvarcount"
    doc: |
      Number of highly variable features to detect. Used for datasets integration,
      scaling, and dimensional reduction.
      Default: 3000

  dimensionality:
    type: int?
    inputBinding:
      prefix: "--ndim"
    doc: |
      Number of principal components to use in UMAP projection and clustering
      (from 1 to 50). Use Elbow plot to adjust this parameter.
      Default: 10

  umap_spread:
    type: float?
    inputBinding:
      prefix: "--spread"
    doc: |
      The effective scale of embedded points on UMAP. In combination with mindist
      this determines how clustered/clumped the embedded points are.
      Default: 1

  umap_mindist:
    type: float?
    inputBinding:
      prefix: "--mindist"
    doc: |
      Controls how tightly the embedding is allowed compress points together on UMAP.
      Larger values ensure embedded points are moreevenly distributed, while smaller
      values allow the algorithm to optimise more accurately with regard to local structure.
      Sensible values are in the range 0.001 to 0.5.
      Default:  0.3

  umap_nneighbors:
    type: int?
    inputBinding:
      prefix: "--nneighbors"
    doc: |
      Determines the number of neighboring points used in UMAP. Larger values will result
      in more global structure being preserved at the loss of detailed local structure.
      In general this parameter should often be in the range 5 to 50.
      Default: 30

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
    inputBinding:
      prefix: "--umetric"
    doc: |
      The metric to use to compute distances in high dimensional space for UMAP.
      Default: cosine

  umap_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "uwot"
      - "uwot-learn"
      - "umap-learn"
    inputBinding:
      prefix: "--umethod"
    doc: |
      UMAP implementation to run.
      Default: uwot

  cluster_metric:
    type:
    - "null"
    - type: enum
      symbols:
      - "euclidean"
      - "cosine"
      - "manhattan"
      - "hamming"
    inputBinding:
      prefix: "--ametric"
    doc: |
      Distance metric used by the nearest neighbors algorithm when running clustering.
      Default: cosine

  resolution:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--resolution"
    doc: |
      Clustering resolution. Can be set as an array.
      Default: 0.4 0.6 0.8 1.0 1.4

  minimum_logfc:
    type: float?
    inputBinding:
      prefix: "--logfc"
    doc: |
      Include only those genes that on average have log fold change difference in
      expression between every tested pair of clusters not lower than this value.
      Default: 0.25

  minimum_pct:
    type: float?
    inputBinding:
      prefix: "--minpct"
    doc: |
      Include only those features that are detected in not lower than this fraction
      of cells in either of the two tested clusters.
      Default: 0.1

  only_positive_markers:
    type: boolean?
    inputBinding:
      prefix: "--onlypos"
    doc: |
      Return only positive markers when running gene markers identification.
      Default: false

  no_sct:
    type: boolean?
    inputBinding:
      prefix: "--nosct"
    doc: |
      Do not use SCTransform when running datasets integration. Use LogNormalize instead.
      Default: false

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
    inputBinding:
      prefix: "--testuse"
    doc: |
      Statistical test to use for gene markers identification.
      Default: wilcox

  species:
    type:
    - "null"
    - type: enum
      symbols:
      - "hs"
      - "mm"
      - "none"
    inputBinding:
      prefix: "--species"
    doc: |
      Select species for gene name conversion when running cell type prediction
      with Garnett classifier.
      Default: do not convert gene names

  export_pdf_plots:
    type: boolean?
    inputBinding:
      prefix: "--pdf"
    doc: |
      Export plots in PDF.
      Default: false

  export_rds_data:
    type: boolean?
    inputBinding:
      prefix: "--rds"
    doc: |
      Save Seurat data to RDS file.
      Default: false

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output prefix.
      Default: ./seurat

  threads:
    type: int?
    inputBinding:
      prefix: "--threads"
    doc: |
      Threads number
      Default: 1


outputs:

  raw_cell_count_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_cell_count.png"
    doc: |
      Number of cells per dataset (not filtered).
      PNG format

  raw_cell_count_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_cell_count.pdf"
    doc: |
      Number of cells per dataset (not filtered).
      PDF format

  raw_umi_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_umi_dnst_spl_by_cond.png"
    doc: |
      Split by condition UMI density per cell (not filtered).
      PNG format

  raw_umi_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_umi_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition UMI density per cell (not filtered).
      PDF format

  raw_gene_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst_spl_by_cond.png"
    doc: |
      Split by condition gene density per cell (not filtered).
      PNG format

  raw_gene_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition gene density per cell (not filtered).
      PDF format

  raw_gene_umi_corr_spl_by_ident_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi_corr_spl_by_ident.png"
    doc: |
      Split by identity genes vs UMIs per cell correlation (not filtered).
      PNG format

  raw_gene_umi_corr_spl_by_ident_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi_corr_spl_by_ident.pdf"
    doc: |
      Split by identity genes vs UMIs per cell correlation (not filtered).
      PDF format

  raw_mito_perc_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_mito_perc_dnst_spl_by_cond.png"
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (not filtered).
      PNG format

  raw_mito_perc_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_mito_perc_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (not filtered).
      PDF format

  raw_nvlt_score_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_score_dnst_spl_by_cond.png"
    doc: |
      Split by condition novelty score density per cell (not filtered).
      PNG format

  raw_nvlt_score_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_score_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition novelty score density per cell (not filtered).
      PDF format

  raw_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_qc_mtrcs.png"
    doc: |
      QC metrics densities per cell (not filtered).
      PNG format

  raw_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_qc_mtrcs.pdf"
    doc: |
      QC metrics densities per cell (not filtered).
      PDF format

  raw_qc_mtrcs_gr_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_qc_mtrcs_gr_by_cond.png"
    doc: |
      Grouped by condition QC metrics densities per cell (not filtered).
      PNG format

  raw_qc_mtrcs_gr_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_qc_mtrcs_gr_by_cond.pdf"
    doc: |
      Grouped by condition QC metrics densities per cell (not filtered).
      PDF format


  fltr_cell_count_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_cell_count.png"
    doc: |
      Number of cells per dataset (filtered).
      PNG format

  fltr_cell_count_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_cell_count.pdf"
    doc: |
      Number of cells per dataset (filtered).
      PDF format

  fltr_umi_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_umi_dnst_spl_by_cond.png"
    doc: |
      Split by condition UMI density per cell (filtered).
      PNG format

  fltr_umi_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_umi_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition UMI density per cell (filtered).
      PDF format

  fltr_gene_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_dnst_spl_by_cond.png"
    doc: |
      Split by condition gene density per cell (filtered).
      PNG format

  fltr_gene_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition gene density per cell (filtered).
      PDF format

  fltr_gene_umi_corr_spl_by_ident_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_umi_corr_spl_by_ident.png"
    doc: |
      Split by identity genes vs UMIs per cell correlation (filtered).
      PNG format

  fltr_gene_umi_corr_spl_by_ident_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_umi_corr_spl_by_ident.pdf"
    doc: |
      Split by identity genes vs UMIs per cell correlation (filtered).
      PDF format

  fltr_mito_perc_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_mito_perc_dnst_spl_by_cond.png"
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (filtered).
      PNG format

  fltr_mito_perc_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_mito_perc_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (filtered).
      PDF format

  fltr_nvlt_score_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_nvlt_score_dnst_spl_by_cond.png"
    doc: |
      Split by condition novelty score density per cell (filtered).
      PNG format

  fltr_nvlt_score_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_nvlt_score_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition novelty score density per cell (filtered).
      PDF format

  fltr_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_qc_mtrcs.png"
    doc: |
      QC metrics densities per cell (filtered).
      PNG format

  fltr_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_qc_mtrcs.pdf"
    doc: |
      QC metrics densities per cell (filtered).
      PDF format

  fltr_qc_mtrcs_gr_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_qc_mtrcs_gr_by_cond.png"
    doc: |
      Grouped by condition QC metrics densities per cell (filtered).
      PNG format

  fltr_qc_mtrcs_gr_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_qc_mtrcs_gr_by_cond.pdf"
    doc: |
      Grouped by condition QC metrics densities per cell (filtered).
      PDF format


  fltr_pca_spl_by_ph_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_pca_spl_by_ph.png"
    doc: |
      Split by cell cycle phase PCA of filtered unintegrated/scaled datasets.
      PNG format

  fltr_pca_spl_by_ph_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_pca_spl_by_ph.pdf"
    doc: |
      Split by cell cycle phase PCA of filtered unintegrated/scaled datasets.
      PDF format

  fltr_pca_spl_by_mito_perc_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_pca_spl_by_mito_perc.png"
    doc: |
      Split by level of transcripts mapped to mitochondrial genes PCA of filtered unintegrated/scaled datasets.
      PNG format

  fltr_pca_spl_by_mito_perc_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_pca_spl_by_mito_perc.pdf"
    doc: |
      Split by level of transcripts mapped to mitochondrial genes PCA of filtered unintegrated/scaled datasets.
      PDF format

  fltr_umap_spl_by_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_umap_spl_by_idnt.png"
    doc: |
      Split by identity UMAP projected PCA of filtered unintegrated/scaled datasets.
      PNG format

  fltr_umap_spl_by_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_umap_spl_by_idnt.pdf"
    doc: |
      Split by identity UMAP projected PCA of filtered unintegrated/scaled datasets.
      PDF format


  ntgr_elbow_plot_png:
    type: File?
    outputBinding:
      glob: "*_ntgr_elbow.png"
    doc: |
      Elbow plot from PCA of filtered integrated/scaled datasets.
      PNG format

  ntgr_elbow_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_ntgr_elbow.pdf"
    doc: |
      Elbow plot from PCA of filtered integrated/scaled datasets.
      PDF format

  ntgr_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_ntgr_pca.png"
    doc: |
      PCA of filtered integrated/scaled datasets.
      PNG format

  ntgr_pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_ntgr_pca.pdf"    
    doc: |
      PCA of filtered integrated/scaled datasets.
      PDF format

  ntgr_pca_heatmap_png:
    type: File?
    outputBinding:
      glob: "*_ntgr_pca_heatmap.png"
    doc: |
      Genes per cells expression heatmap sorted by their PC scores from PCA of filtered integrated/scaled datasets.
      PNG format

  ntgr_pca_heatmap_pdf:
    type: File?
    outputBinding:
      glob: "*_ntgr_pca_heatmap.pdf"
    doc: |
      Genes per cells expression heatmap sorted by their PC scores from PCA of filtered integrated/scaled datasets.
      PDF format

  ntgr_pca_loadings_plot_png:
    type: File?
    outputBinding:
      glob: "*_ntgr_pca_loadings.png"
    doc: |
      PC scores of the most variant genes from PCA of filtered integrated/scaled datasets.
      PNG format

  ntgr_pca_loadings_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_ntgr_pca_loadings.pdf"
    doc: |
      PC scores of the most variant genes from PCA of filtered integrated/scaled datasets.
      PDF format

  ntgr_umap_spl_by_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_ntgr_umap_spl_by_idnt.png"
    doc: |
      Split by identity UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format

  ntgr_umap_spl_by_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_ntgr_umap_spl_by_idnt.pdf"
    doc: |
      Split by identity UMAP projected PCA of filtered integrated/scaled datasets.
      PDF format


  clst_umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_res_*.png"
    doc: |
      Clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format

  clst_umap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_res_*.pdf"
    doc: |
      Clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PDF format

  clst_umap_spl_by_cond_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_spl_by_cond_res_*.png"
    doc: |
      Split by condition clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format

  clst_umap_spl_by_cond_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_spl_by_cond_res_*.pdf"
    doc: |
      Split by condition clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PDF format

  clst_umap_ctype_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_ctype_res_*.png"
    doc: |
      Grouped by predicted cell types UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format

  clst_umap_ctype_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_ctype_res_*.pdf"
    doc: |
      Grouped by predicted cell types UMAP projected PCA of filtered integrated/scaled datasets.
      PDF format

  clst_umap_spl_by_ph_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_spl_by_ph_res_*.png"
    doc: |
      Split by cell cycle phase clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format

  clst_umap_spl_by_ph_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_spl_by_ph_res_*.pdf"
    doc: |
      Split by cell cycle phase clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PDF format

  clst_qc_mtrcs_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_qc_mtrcs_res_*.png"
    doc: |
      QC metrics for clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format

  clst_qc_mtrcs_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_qc_mtrcs_res_*.pdf"
    doc: |
      QC metrics for clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PDF format

  expr_avg_per_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_avg_per_clst_res_*.png"
    doc: |
      Scaled average log normalized gene expression per cluster of filtered integrated/scaled datasets.
      PNG format

  expr_avg_per_clst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_avg_per_clst_res_*.pdf"
    doc: |
      Scaled average log normalized gene expression per cluster of filtered integrated/scaled datasets.
      PDF format

  expr_per_clst_cell_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_per_clst_cell_res_*.png"
    doc: |
      Log normalized gene expression per cell of clustered filtered integrated/scaled datasets.
      PNG format

  expr_per_clst_cell_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_per_clst_cell_res_*.pdf"
    doc: |
      Log normalized gene expression per cell of clustered filtered integrated/scaled datasets.
      PDF format

  expr_clst_heatmap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_clst_heatmap_res_*.png"
    doc: |
      Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets.
      PNG format

  expr_clst_heatmap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_clst_heatmap_res_*.pdf"
    doc: |
      Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets.
      PDF format

  expr_dnst_per_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_dnst_per_clst_res_*.png"
    doc: |
      Log normalized gene expression densities per cluster of filtered integrated/scaled datasets.
      PNG format

  expr_dnst_per_clst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_dnst_per_clst_res_*.pdf"
    doc: |
      Log normalized gene expression densities per cluster of filtered integrated/scaled datasets.
      PDF format

  expr_avg_per_ctype_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_avg_per_ctype_res_*.png"
    doc: |
      Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets.
      PNG format

  expr_avg_per_ctype_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_avg_per_ctype_res_*.pdf"
    doc: |
      Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets.
      PDF format

  expr_per_ctype_cell_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_per_ctype_cell_res_*.png"
    doc: |
      Log normalized gene expression per cell of clustered filtered integrated/scaled datasets with predicted cell types.
      PNG format

  expr_per_ctype_cell_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_per_ctype_cell_res_*.pdf"
    doc: |
      Log normalized gene expression per cell of clustered filtered integrated/scaled datasets with predicted cell types.
      PDF format

  expr_ctype_heatmap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_ctype_heatmap_res_*.png"
    doc: |
      Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets with predicted cell types.
      PNG format

  expr_ctype_heatmap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_ctype_heatmap_res_*.pdf"
    doc: |
      Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets with predicted cell types.
      PDF format

  expr_dnst_per_ctype_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_dnst_per_ctype_res_*.png"
    doc: |
      Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets.
      PNG format

  expr_dnst_per_ctype_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_dnst_per_ctype_res_*.pdf"
    doc: |
      Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets.
      PDF format

  clst_pttv_gene_markers:
    type: File
    outputBinding:
      glob: "*_clst_pttv_gene_markers.tsv"
    doc: |
      Putative gene markers file for all clusters and all resolutions.
      TSV format

  clst_csrvd_gene_markers:
    type: File
    outputBinding:
      glob: "*_clst_csrvd_gene_markers.tsv"
    doc: |
      Conserved gene markers file for all clusters and all resolutions.
      TSV format

  seurat_clst_data_rds:
    type: File?
    outputBinding:
      glob: "*_clst_data.rds"
    doc: |
      Clustered filtered integrated/scaled Seurat data.
      RDS format

  cellbrowser_config_data:
    type: Directory
    outputBinding:
      glob: "*_cellbrowser"
    doc: |
      Directory with UCSC Cellbrowser configuration data

  cellbrowser_html_data:
    type: Directory
    outputBinding:
      glob: "*_cellbrowser/html_data"
    doc: |
      Directory with UCSC Cellbrowser formatted html data

  cellbrowser_html_file:
    type: File
    outputBinding:
      glob: "*_cellbrowser/html_data/index.html"
    doc: |
      HTML index file from the directory with UCSC Cellbrowser formatted html data


  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["run_seurat.R"]


stdout: seurat_cluster_stdout.log
stderr: seurat_cluster_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Seurat cluster"
s:name: "Seurat cluster"
s:alternateName: "Runs Seurat for comparative scRNA-seq analysis of across experimental conditions"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/seurat-cluster.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
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
  Seurat cluster
  ==============

  The joint analysis of multiple scRNA-Seq datasets with [Seurat](https://satijalab.org/seurat/) starts with evaluation of common
  single-cell quality control (QC) metrics – genes and UMIs counts, percentage of mitochondrial genes
  expressed. QC allows to get a general overview of the datasets quality as well as to define filtering
  thresholds for dead or low-quality cells removal. Filtered merged datasets are then being processed
  with the integration algorithm. Its main goal is to identify integration anchors – pairs of cells that can
  “pull together” the same cell type populations from the different datasets. An integration algorithm
  can also solve batch correction problem by regressing out the unwanted sources of variation. The
  integrated data then undergo the dimensionality reduction processing that starts from the principal
  component analysis (PCA). Based on the PCA results the uniform manifold approximation and
  projection (UMAP) and clustering analysis are run with the principal components of the highest
  variance. Clustered data are then used for gene markers identification. These genes are differentially
  expressed between clusters and can be used for cell types assignment.
  More details about scRNA-Seq integration analysis with Seurat can be found in the official
  [documentation](https://satijalab.org/seurat/articles/integration_introduction.html).


s:about: |
  usage: run_seurat.R [-h] --mex MEX [MEX ...] --identity
                                    IDENTITY [--condition CONDITION]
                                    [--classifier CLASSIFIER]
                                    [--cellcycle CELLCYCLE]
                                    [--barcodes BARCODES] [--mincells MINCELLS]
                                    [--minfeatures [MINFEATURES [MINFEATURES ...]]]
                                    [--maxfeatures [MAXFEATURES [MAXFEATURES ...]]]
                                    [--minumi [MINUMI [MINUMI ...]]]
                                    [--minnovelty [MINNOVELTY [MINNOVELTY ...]]]
                                    [--maxmt MAXMT] [--mitopattern MITOPATTERN]
                                    [--features [FEATURES [FEATURES ...]]]
                                    [--regresscellcycle] [--regressmt]
                                    [--highvarcount HIGHVARCOUNT] [--ndim NDIM]
                                    [--resolution [RESOLUTION [RESOLUTION ...]]]
                                    [--logfc LOGFC] [--minpct MINPCT]
                                    [--onlypos]
                                    [--testuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}]
                                    [--species {hs,mm,none}] [--pdf] [--rds]
                                    [--output OUTPUT] [--threads THREADS]

  Runs Seurat for comparative scRNA-seq analysis of across experimental
  conditions

  optional arguments:
    -h, --help            show this help message and exit
    --mex MEX [MEX ...]   Path to the folder with not normalized aggregated
                          feature-barcode matrix from Cell Ranger Aggregate in
                          MEX format. If multiple locations provided data is
                          assumed to be not aggregated (outputs from multiple
                          Cell Ranger Count runs) and will be merged.
    --identity IDENTITY   Path to the metadata TSV/CSV file to set the datasets
                          identities. If --mex points to the Cell Ranger
                          Aggregate outputs, the aggregation.csv file can be
                          used as well. If multiple locations were provided
                          through --mex, the file should include at least one
                          column - 'library_id', and be sorted based on the the
                          order of locations provided in --mex.
    --condition CONDITION
                          Path to the TSV/CSV file to define datasets grouping.
                          First column - 'library_id' with the values provided
                          in the correspondent column of the --identity file,
                          second column 'condition'. Default: each dataset is
                          assigned to a separate group.
    --classifier CLASSIFIER
                          Path to the Garnett classifier RDS file for cell type
                          prediction. Default: skip cell type prediction.
    --cellcycle CELLCYCLE
                          Path to the TSV/CSV file with cell cycle data. First
                          column - 'phase', second column 'gene_id'. Default:
                          skip cell cycle score assignment.
    --barcodes BARCODES   Path to the headerless TSV/CSV file with the list of
                          barcodes to select cells of interest (one barcode per
                          line). Prefilters input feature-barcode matrix to
                          include only selected cells. Default: use all cells.
    --mincells MINCELLS   Include only features detected in at least this many
                          cells. Applied to aggregated feature-barcode matrix
                          from Cell Ranger Aggregate. Ignored when --mex points
                          to the locations of multiple Cell Ranger Count runs.
                          Default: 5
    --minfeatures [MINFEATURES [MINFEATURES ...]]
                          Include cells where at least this many features are
                          detected. If multiple values provided each of them
                          will be applied to the correspondent dataset from the
                          --mex input. Default: 250 (applied to all datasets)
    --maxfeatures [MAXFEATURES [MAXFEATURES ...]]
                          Include cells with the number of features not bigger
                          than this value. If multiple values provided each of
                          them will be applied to the correspondent dataset from
                          the --mex input. Default: 5000 (applied to all
                          datasets)
    --minumi [MINUMI [MINUMI ...]]
                          Include cells where at least this many UMIs
                          (transcripts) are detected. If multiple values
                          provided each of them will be applied to the
                          correspondent dataset from the --mex input. Default:
                          500 (applied to all datasets)
    --minnovelty [MINNOVELTY [MINNOVELTY ...]]
                          Include cells with the novelty score not lower than
                          this value, calculated as log10(genes)/log10(UMIs). If
                          multiple values provided each of them will be applied
                          to the correspondent dataset from the --mex input.
                          Default: 0.8 (applied to all datasets)
    --maxmt MAXMT         Include cells with the percentage of transcripts
                          mapped to mitochondrial genes not bigger than this
                          value. Default: 5
    --mitopattern MITOPATTERN
                          Regex pattern to identify mitochondrial genes.
                          Default: '^Mt-'
    --features [FEATURES [FEATURES ...]]
                          Features of interest to evaluate expression. Default:
                          None
    --regresscellcycle    Regress cell cycle as a confounding source of
                          variation. Default: false
    --regressmt           Regress mitochondrial genes expression as a
                          confounding source of variation. Default: false
    --highvarcount HIGHVARCOUNT
                          Number of highly variable features to detect. Used for
                          datasets integration, scaling, and dimensional
                          reduction. Default: 3000
    --ndim NDIM           Number of principal components to use in UMAP
                          projection and clustering (from 1 to 50). Use Elbow
                          plot to adjust this parameter. Default: 10
    --resolution [RESOLUTION [RESOLUTION ...]]
                          Clustering resolution. Can be set as an array.
                          Default: 0.4 0.6 0.8 1.0 1.4
    --logfc LOGFC         Include only those genes that on average have log fold
                          change difference in expression between every tested
                          pair of clusters not lower than this value. Default:
                          0.25
    --minpct MINPCT       Include only those features that are detected in not
                          lower than this fraction of cells in either of the two
                          tested clusters. Default: 0.1
    --onlypos             Return only positive markers when running gene markers
                          identification. Default: false
    --testuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}
                          Statistical test to use for gene markers
                          identification. Default: wilcox
    --species {hs,mm,none}
                          Select species for gene name conversion when running
                          cell type prediction with Garnett classifier. Default:
                          do not convert gene names
    --pdf                 Export plots in PDF. Default: false
    --rds                 Save Seurat data to RDS file. Default: false
    --output OUTPUT       Output prefix. Default: ./seurat
    --threads THREADS     Threads. Default: 1
