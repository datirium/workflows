cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  sc_rnaseq_aggr_sample:
  - "cellranger-aggr.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  filtered_feature_bc_matrix_folder:
    type: File
    label: "scRNA-Seq Cellranger Aggregate Experiment"
    doc: |
      Compressed folder with aggregated filtered feature-barcode matrices in MEX format
    'sd:upstreamSource': "sc_rnaseq_aggr_sample/filtered_feature_bc_matrix_folder"
    'sd:localLabel': true

  aggregation_metadata:
    type: File
    label: "scRNA-Seq Cellranger Aggregate Experiment"
    doc: |
      Aggregation metadata in CSV format
    'sd:upstreamSource': "sc_rnaseq_aggr_sample/aggregation_metadata"
    'sd:localLabel': true

  cell_cycle_data:
    type: File
    label: "TSV/CSV file with cell cycle data with 'phase' and 'gene_id' columns"
    doc: |
      TSV/CSV file with cell cycle data. First column - 'phase', second column 'gene_id'

  conditions_data:
    type: File
    label: "TSV/CSV file to define datasets conditions with 'library_id' and 'condition' columns"
    doc: |
      Path to the TSV/CSV file to define datasets conditions
      for grouping. First column - 'library_id' with values
      from the --identity file, second column 'condition'.
      Default: each dataset is assigned to its own biological
      condition

  classifier_rds:
    type: File
    label: "Garnett classifier rds file for cell type prediction"
    doc: |
      Path to the Garnett classifier rds file for cell type prediction

  species:
    type:
      type: enum
      symbols:
      - "hs"
      - "mm"
    label: "Species for gene name conversion"
    doc: |
      Select species for gene name conversion when running cell
      type prediction. Either "hs" or "mm"

  minimum_cells:
    type: int?
    default: 10
    label: "Include features detected in at least this many cells"
    doc: |
      Include features detected in at least this many cells
      (applied to thoughout all datasets together).
    'sd:layout':
      advanced: true

  minimum_features:
    type: int?
    default: 250
    label: "Include cells where at least this many features are detected"
    doc: |
      Include cells where at least this many features are detected.
    'sd:layout':
      advanced: true

  minimum_umis:
    type: int?
    default: 500
    label: "Include cells where at least this many UMI are detected"
    doc: |
      Include cells where at least this many UMI are detected.
    'sd:layout':
      advanced: true

  minimum_novelty_score:
    type: float?
    default: 0.8
    label: "Include cells with the novelty score not lower that this value"
    doc: |
      Include cells with the novelty score not lower that this
      value (calculated as log10(genes)/log10(UMIs)).
    'sd:layout':
      advanced: true

  maximum_mito_perc:
    type: float?
    default: 5
    label: "Include cells with the mitochondrial contamination percentage not bigger that this value"
    doc: |
      Include cells with the mitochondrial contamination percentage
      not bigger that this value.
    'sd:layout':
      advanced: true

  mito_pattern:
    type: string?
    default: "^Mt-"
    label: "Pattern to identify mitochondrial reads"
    doc: |
      Regex pattern to identify mitochondrial reads.
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

  high_var_features_count:
    type: int?
    default: 3000
    label: "Number of higly variable features to detect"
    doc: |
      Number of higly variable features to detect.
    'sd:layout':
      advanced: true

  dimensionality:
    type: int?
    default: 10
    label: "Number of principal components to use in clustering (1:50)"
    doc: |
      Number of principal components to use in clustering (1:50).
      Use Elbow plot to adjust this parameter.
    'sd:layout':
      advanced: true

  resolution:
    type: float?
    default: 0.1
    label: "Clustering resolution"
    doc: |
      Clustering resolution
    'sd:layout':
      advanced: true

  minimum_logfc:
    type: float?
    default: 0.25
    label: "Log fold change threshold for conserved gene markers identification"
    doc: |
      Log fold change threshold for conserved gene markers identification.
    'sd:layout':
      advanced: true

  minimum_pct:
    type: float?
    default: 0.1
    label: "Minimum fraction of cells where genes used for gene markers identification should be detected"
    doc: |
      Minimum fraction of cells where genes used for conserved gene markers
      identification should be detected in either of two tested clusters.
    'sd:layout':
      advanced: true

  only_positive_markers:
    type: boolean?
    default: false
    label: "Return only positive markers when running conserved gene markers identification"
    doc: |
      Return only positive markers when running conserved gene markers
      identification.
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


outputs:

  raw_cell_count_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_cell_count_plot_png
    label: "Raw number of cells per dataset plot"
    doc: |
      Raw number of cells per dataset plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Raw QC'
        Caption: 'Raw number of cells per dataset'

  raw_cell_count_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_cell_count_plot_pdf
    label: "Raw number of cells per dataset plot"
    doc: |
      Raw number of cells per dataset plot in PDF format

  raw_umi_log10_density_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_umi_log10_density_plot_png
    label: "Raw UMI per cell log10 density plot"
    doc: |
      Raw UMI per cell log10 density plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Raw QC'
        Caption: 'Raw UMI per cell log10 density'

  raw_umi_log10_density_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_umi_log10_density_plot_pdf
    label: "Raw UMI per cell log10 density plot"
    doc: |
      Raw UMI per cell log10 density plot in PDF format

  raw_gene_log10_density_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_gene_log10_density_plot_png
    label: "Raw gene per cell log10 density plot"
    doc: |
      Raw gene per cell log10 density plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Raw QC'
        Caption: 'Raw gene per cell log10 density'

  raw_gene_log10_density_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_gene_log10_density_plot_pdf
    label: "Raw gene per cell log10 density plot"
    doc: |
      Raw gene per cell log10 density plot in PDF format

  raw_gene_umi_log10_correlation_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_gene_umi_log10_correlation_plot_png
    label: "Raw log10 gene vs log10 UMI correlation plot"
    doc: |
      Raw log10 gene vs log10 UMI correlation plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Raw QC'
        Caption: 'Raw log10 gene vs log10 UMI correlation'

  raw_gene_umi_log10_correlation_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_gene_umi_log10_correlation_plot_pdf
    label: "Raw log10 gene vs log10 UMI correlation plot"
    doc: |
      Raw log10 gene vs log10 UMI correlation plot in PDF format

  raw_mito_perc_log10_density_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_mito_perc_log10_density_plot_png
    label: "Raw mitochondrial gene percentage per cell log10 density plot"
    doc: |
      Raw mitochondrial gene percentage per cell log10 density plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Raw QC'
        Caption: 'Raw mitochondrial gene percentage per cell log10 density'

  raw_mito_perc_log10_density_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_mito_perc_log10_density_plot_pdf
    label: "Raw mitochondrial gene percentage per cell log10 density plot"
    doc: |
      Raw mitochondrial gene percentage per cell log10 density plot in PDF format

  raw_novelty_score_log10_density_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_novelty_score_log10_density_plot_png
    label: "Raw novelty score (log10Gene/log10UMI) per cell log10 density plot"
    doc: |
      Raw novelty score (log10Gene/log10UMI) per cell log10 density plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Raw QC'
        Caption: 'Raw novelty score (log10Gene/log10UMI) per cell log10 density'

  raw_novelty_score_log10_density_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_novelty_score_log10_density_plot_pdf
    label: "Raw novelty score (log10Gene/log10UMI) per cell log10 density plot"
    doc: |
      Raw novelty score (log10Gene/log10UMI) per cell log10 density plot in PDF format

  raw_qc_metrics_vln_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_qc_metrics_vln_plot_png
    label: "Raw QC metrics violin plot"
    doc: |
      Raw QC metrics violin plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Raw QC'
        Caption: 'Raw QC metrics violin plot'

  raw_qc_metrics_vln_plot_pdf:
    type: File?
    outputSource: seurat_cluster/raw_qc_metrics_vln_plot_pdf
    label: "Raw QC metrics violin plot"
    doc: |
      Raw QC metrics violin plot in PDF format

  raw_qc_metrics_vln_plot_gr_by_cond_png:
    type: File?
    outputSource: seurat_cluster/raw_qc_metrics_vln_plot_gr_by_cond_png
    label: "Raw QC metrics violin plot grouped by condition"
    doc: |
      Raw QC metrics violin plot grouped by condition in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Raw QC'
        Caption: 'Raw QC metrics violin plot grouped by condition'

  raw_qc_metrics_vln_plot_gr_by_cond_pdf:
    type: File?
    outputSource: seurat_cluster/raw_qc_metrics_vln_plot_gr_by_cond_pdf
    label: "Raw QC metrics violin plot grouped by condition"
    doc: |
      Raw QC metrics violin plot grouped by condition in PDF format

  filt_cell_count_plot_png:
    type: File?
    outputSource: seurat_cluster/filt_cell_count_plot_png
    label: "Filtered number of cells per dataset plot"
    doc: |
      Filtered number of cells per dataset plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Filtered number of cells per dataset'

  filt_cell_count_plot_pdf:
    type: File?
    outputSource: seurat_cluster/filt_cell_count_plot_pdf
    label: "Filtered number of cells per dataset plot"
    doc: |
      Filtered number of cells per dataset plot in PDF format

  filt_umi_log10_density_plot_png:
    type: File?
    outputSource: seurat_cluster/filt_umi_log10_density_plot_png
    label: "Filtered UMI per cell log10 density plot"
    doc: |
      Filtered UMI per cell log10 density plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Filtered UMI per cell log10 density'

  filt_umi_log10_density_plot_pdf:
    type: File?
    outputSource: seurat_cluster/filt_umi_log10_density_plot_pdf
    label: "Filtered UMI per cell log10 density plot"
    doc: |
      Filtered UMI per cell log10 density plot in PDF format

  filt_gene_log10_density_plot_png:
    type: File?
    outputSource: seurat_cluster/filt_gene_log10_density_plot_png
    label: "Filtered gene per cell log10 density plot"
    doc: |
      Filtered gene per cell log10 density plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Filtered gene per cell log10 density'

  filt_gene_log10_density_plot_pdf:
    type: File?
    outputSource: seurat_cluster/filt_gene_log10_density_plot_pdf
    label: "Filtered gene per cell log10 density plot"
    doc: |
      Filtered gene per cell log10 density plot in PDF format

  filt_gene_umi_log10_correlation_plot_png:
    type: File?
    outputSource: seurat_cluster/filt_gene_umi_log10_correlation_plot_png
    label: "Filtered log10 gene vs log10 UMI correlation plot"
    doc: |
      Filtered log10 gene vs log10 UMI correlation plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Filtered log10 gene vs log10 UMI correlation'

  filt_gene_umi_log10_correlation_plot_pdf:
    type: File?
    outputSource: seurat_cluster/filt_gene_umi_log10_correlation_plot_pdf
    label: "Filtered log10 gene vs log10 UMI correlation plot"
    doc: |
      Filtered log10 gene vs log10 UMI correlation plot in PDF format

  filt_mito_perc_log10_density_plot_png:
    type: File?
    outputSource: seurat_cluster/filt_mito_perc_log10_density_plot_png
    label: "Filtered mitochondrial gene percentage per cell log10 density plot"
    doc: |
      Filtered mitochondrial gene percentage per cell log10 density plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Filtered mitochondrial gene percentage per cell log10 density'

  filt_mito_perc_log10_density_plot_pdf:
    type: File?
    outputSource: seurat_cluster/filt_mito_perc_log10_density_plot_pdf
    label: "Filtered mitochondrial gene percentage per cell log10 density plot"
    doc: |
      Filtered mitochondrial gene percentage per cell log10 density plot in PDF format

  filt_novelty_score_log10_density_plot_png:
    type: File?
    outputSource: seurat_cluster/filt_novelty_score_log10_density_plot_png
    label: "Filtered novelty score (log10Gene/log10UMI) per cell log10 density plot"
    doc: |
      Filtered novelty score (log10Gene/log10UMI) per cell log10 density plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Filtered novelty score (log10Gene/log10UMI) per cell log10 density'

  filt_novelty_score_log10_density_plot_pdf:
    type: File?
    outputSource: seurat_cluster/filt_novelty_score_log10_density_plot_pdf
    label: "Filtered novelty score (log10Gene/log10UMI) per cell log10 density plot"
    doc: |
      Filtered novelty score (log10Gene/log10UMI) per cell log10 density plot in PDF format

  filt_qc_metrics_vln_plot_png:
    type: File?
    outputSource: seurat_cluster/filt_qc_metrics_vln_plot_png
    label: "Filtered QC metrics violin plot"
    doc: |
      Filtered QC metrics violin plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Filtered QC metrics violin'

  filt_qc_metrics_vln_plot_pdf:
    type: File?
    outputSource: seurat_cluster/filt_qc_metrics_vln_plot_pdf
    label: "Filtered QC metrics violin plot"
    doc: |
      Filtered QC metrics violin plot in PDF format

  filt_qc_metrics_vln_plot_gr_by_cond_png:
    type: File?
    outputSource: seurat_cluster/filt_qc_metrics_vln_plot_gr_by_cond_png
    label: "Filtered QC metrics violin plot grouped by condition"
    doc: |
      Filtered QC metrics violin plot grouped by condition in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Filtered QC metrics violin plot grouped by condition'

  filt_qc_metrics_vln_plot_gr_by_cond_pdf:
    type: File?
    outputSource: seurat_cluster/filt_qc_metrics_vln_plot_gr_by_cond_pdf
    label: "Filtered QC metrics violin plot grouped by condition"
    doc: |
      Filtered QC metrics violin plot grouped by condition in PDF format

  filt_unint_cell_cycle_eff_pca_plot_png:
    type: File?
    outputSource: seurat_cluster/filt_unint_cell_cycle_eff_pca_plot_png
    label: "Filtered unintegrated PCA plot to evaluate cell cycle as a source of unwanted variation"
    doc: |
      Filtered unintegrated PCA plot to evaluate cell cycle as a source of
      unwanted variation in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Unwanted Variation'
        Caption: 'Filtered unintegrated PCA plot to evaluate cell cycle as a source of unwanted variation'

  filt_unint_cell_cycle_eff_pca_plot_pdf:
    type: File?
    outputSource: seurat_cluster/filt_unint_cell_cycle_eff_pca_plot_pdf
    label: "Filtered unintegrated PCA plot to evaluate cell cycle as a source of unwanted variation"
    doc: |
      Filtered unintegrated PCA plot to evaluate cell cycle as a source of
      unwanted variation in PDF format

  filt_unint_mito_perc_eff_pca_plot_png:
    type: File?
    outputSource: seurat_cluster/filt_unint_mito_perc_eff_pca_plot_png
    label: "Filtered unintegrated PCA plot to evaluate mitochondrial contamination percentage as a source of unwanted variation"
    doc: |
      Filtered unintegrated PCA plot to evaluate mitochondrial contamination
      percentage as a source of unwanted variation in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Unwanted Variation'
        Caption: 'Filtered unintegrated PCA plot to evaluate mitochondrial contamination percentage as a source of unwanted variation'

  filt_unint_mito_perc_eff_pca_plot_pdf:
    type: File?
    outputSource: seurat_cluster/filt_unint_mito_perc_eff_pca_plot_pdf
    label: "Filtered unintegrated PCA plot to evaluate mitochondrial contamination percentage as a source of unwanted variation"
    doc: |
      Filtered unintegrated PCA plot to evaluate mitochondrial contamination
      percentage as a source of unwanted variation in PDF format

  filt_unint_umap_plot_spl_by_ident_png:
    type: File?
    outputSource: seurat_cluster/filt_unint_umap_plot_spl_by_ident_png
    label: "Filtered unintegrated UMAP plot split by cell identity"
    doc: |
      Filtered unintegrated UMAP plot split by cell identity in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Unwanted Variation'
        Caption: 'Filtered unintegrated UMAP plot split by cell identity'

  filt_unint_umap_plot_spl_by_ident_pdf:
    type: File?
    outputSource: seurat_cluster/filt_unint_umap_plot_spl_by_ident_pdf
    label: "Filtered unintegrated UMAP plot split by cell identity"
    doc: |
      Filtered unintegrated UMAP plot split by cell identity in PDF format

  filt_int_elbow_plot_png:
    type: File?
    outputSource: seurat_cluster/filt_int_elbow_plot_png
    label: "Filtered integrated Elbow plot to evaluate data dimensionality"
    doc: |
      Filtered integrated Elbow plot to evaluate data dimensionality
      in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Dimensionality Evaluation'
        Caption: 'Filtered integrated Elbow plot to evaluate data dimensionality'

  filt_int_elbow_plot_pdf:
    type: File?
    outputSource: seurat_cluster/filt_int_elbow_plot_pdf
    label: "Filtered integrated Elbow plot to evaluate data dimensionality"
    doc: |
      Filtered integrated Elbow plot to evaluate data dimensionality
      in PDF format

  filt_int_pca_plot_png:
    type: File?
    outputSource: seurat_cluster/filt_int_pca_plot_png
    label: "Filtered integrated PCA plot"
    doc: |
      Filtered integrated PCA plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Dimensionality Evaluation'
        Caption: 'Filtered integrated PCA'

  filt_int_pca_plot_pdf:
    type: File?
    outputSource: seurat_cluster/filt_int_pca_plot_pdf
    label: "Filtered integrated PCA plot"
    doc: |
      Filtered integrated PCA plot in PDF format

  filt_int_pca_heatmap_png:
    type: File?
    outputSource: seurat_cluster/filt_int_pca_heatmap_png
    label: "Filtered integrated PCA heatmap to evaluate data dimensionality"
    doc: |
      Filtered integrated PCA heatmap to evaluate data dimensionality
      in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Dimensionality Evaluation'
        Caption: 'Filtered integrated PCA heatmap to evaluate data dimensionality'

  filt_int_pca_heatmap_pdf:
    type: File?
    outputSource: seurat_cluster/filt_int_pca_heatmap_pdf
    label: "Filtered integrated PCA heatmap to evaluate data dimensionality"
    doc: |
      Filtered integrated PCA heatmap to evaluate data dimensionality
      in PDF format

  filt_int_pca_loadings_plot_png:
    type: File?
    outputSource: seurat_cluster/filt_int_pca_loadings_plot_png
    label: "Filtered integrated PCA loadings plot to evaluate data dimensionality"
    doc: |
      Filtered integrated PCA loadings plot to evaluate data dimensionality
      in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Dimensionality Evaluation'
        Caption: 'Filtered integrated PCA loadings plot to evaluate data dimensionality'

  filt_int_pca_loadings_plot_pdf:
    type: File?
    outputSource: seurat_cluster/filt_int_pca_loadings_plot_pdf
    label: "Filtered integrated PCA loadings plot to evaluate data dimensionality"
    doc: |
      Filtered integrated PCA loadings plot to evaluate data dimensionality
      in PDF format

  filt_int_umap_plot_spl_by_ident_png:
    type: File?
    outputSource: seurat_cluster/filt_int_umap_plot_spl_by_ident_png
    label: "Filtered integrated UMAP plot"
    doc: |
      Filtered integrated UMAP plot in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Integration'
        Caption: 'Filtered integrated UMAP'

  filt_int_umap_plot_spl_by_ident_pdf:
    type: File?
    outputSource: seurat_cluster/filt_int_umap_plot_spl_by_ident_pdf
    label: "Filtered integrated UMAP plot"
    doc: |
      Filtered integrated UMAP plot in PDF format

  filt_int_cl_umap_plot_spl_by_cond_res_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/filt_int_cl_umap_plot_spl_by_cond_res_png
    label: "Filtered integrated clustered UMAP plots with variable resolution"
    doc: |
      Filtered integrated clustered UMAP plots with variable resolution
      in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Clustering'
        Caption: 'Filtered integrated clustered UMAP plots with variable resolution'

  filt_int_cl_umap_plot_spl_by_cond_res_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/filt_int_cl_umap_plot_spl_by_cond_res_pdf
    label: "Filtered integrated clustered UMAP plots with variable resolution"
    doc: |
      Filtered integrated clustered UMAP plots with variable resolution
      in PDF format

  filt_int_cl_umap_ctype_pred_plot_res_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/filt_int_cl_umap_ctype_pred_plot_res_png
    label: "Filtered integrated clustered cell type prediction UMAP plots with variable resolution"
    doc: |
      Filtered integrated clustered cell type prediction UMAP plots with
      variable resolution in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Clustering'
        Caption: 'Filtered integrated clustered cell type prediction UMAP plots with variable resolution'

  filt_int_cl_umap_ctype_pred_plot_res_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/filt_int_cl_umap_ctype_pred_plot_res_pdf
    label: "Filtered integrated clustered cell type prediction UMAP plots with variable resolution"
    doc: |
      Filtered integrated clustered cell type prediction UMAP plots with
      variable resolution in PDF format

  filt_int_cl_umap_plot_spl_by_ph_res_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/filt_int_cl_umap_plot_spl_by_ph_res_png
    label: "Filtered integrated clustered UMAP plots split by cell cycle phase with variable resolution"
    doc: |
      Filtered integrated clustered UMAP plots split by cell cycle phase
      with variable resolution in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Post Integration QC'
        Caption: 'Filtered integrated clustered UMAP plots split by cell cycle phase with variable resolution'

  filt_int_cl_umap_plot_spl_by_ph_res_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/filt_int_cl_umap_plot_spl_by_ph_res_pdf
    label: "Filtered integrated clustered UMAP plots split by cell cycle phase with variable resolution"
    doc: |
      Filtered integrated clustered UMAP plots split by cell cycle phase
      with variable resolution in PDF format

  filt_int_cl_umap_qc_metrics_plot_res_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/filt_int_cl_umap_qc_metrics_plot_res_png
    label: "Filtered integrated clustered UMAP QC metrics plots with variable resolution"
    doc: |
      Filtered integrated clustered UMAP QC metrics plots with variable
      resolution in PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Post Integration QC'
        Caption: 'Filtered integrated clustered UMAP QC metrics plots with variable resolution'

  filt_int_cl_umap_qc_metrics_plot_res_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/filt_int_cl_umap_qc_metrics_plot_res_pdf
    label: "Filtered integrated clustered UMAP QC metrics plots with variable resolution"
    doc: |
      Filtered integrated clustered UMAP QC metrics plots with variable
      resolution in PDF format

  seurat_data_rds:
    type: File
    outputSource: seurat_cluster/seurat_data_rds
    label: "Filtered integrated clustered Seurat data in RDS format"
    doc: |
      Filtered integrated clustered Seurat data in RDS format

  conserved_gene_markers:
    type: File
    outputSource: seurat_cluster/conserved_gene_markers
    label: "Conserved gene markers file for all clusters and all resolutions irrespective of condition in TSV format"
    doc: |
      Conserved gene markers file for all clusters and all resolutions
      irrespective of condition in TSV format
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Gene Markers'
        Title: 'Conserved gene markers irrespective of condition'

  compressed_cellbrowser_config_data:
    type: File
    outputSource: compress_cellbrowser_config_data/compressed_folder
    label: "Compressed directory with UCSC Cellbrowser configuration data"
    doc: |
      Directory with UCSC Cellbrowser configuration data

  cellbrowser_html_data:
    type: Directory
    outputSource: seurat_cluster/cellbrowser_html_data
    label: "Directory with UCSC Cellbrowser formatted html data"
    doc: |
      Directory with UCSC Cellbrowser formatted html data

  cellbrowser_html_file:
    type: File
    outputSource: seurat_cluster/cellbrowser_html_file
    label: "HTML index file from the directory with UCSC Cellbrowser formatted html data"
    doc: |
      HTML index file from the directory with UCSC Cellbrowser
      formatted html data
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
      minimum_cells: minimum_cells
      minimum_features: minimum_features
      minimum_umis: minimum_umis
      minimum_novelty_score: minimum_novelty_score
      maximum_mito_perc: maximum_mito_perc
      mito_pattern: mito_pattern
      regress_cellcycle: regress_cellcycle
      regress_mito_perc: regress_mito_perc
      high_var_features_count: high_var_features_count
      dimensionality: dimensionality
      resolution: resolution
      minimum_logfc: minimum_logfc
      minimum_pct: minimum_pct
      only_positive_markers: only_positive_markers
      export_pdf_plots:
        default: true
      export_rds_data:
        default: true
      threads: threads
    out:
    - raw_cell_count_plot_png
    - raw_cell_count_plot_pdf
    - raw_umi_log10_density_plot_png
    - raw_umi_log10_density_plot_pdf
    - raw_gene_log10_density_plot_png
    - raw_gene_log10_density_plot_pdf
    - raw_gene_umi_log10_correlation_plot_png
    - raw_gene_umi_log10_correlation_plot_pdf
    - raw_mito_perc_log10_density_plot_png
    - raw_mito_perc_log10_density_plot_pdf
    - raw_novelty_score_log10_density_plot_png
    - raw_novelty_score_log10_density_plot_pdf
    - raw_qc_metrics_vln_plot_png
    - raw_qc_metrics_vln_plot_pdf
    - raw_qc_metrics_vln_plot_gr_by_cond_png
    - raw_qc_metrics_vln_plot_gr_by_cond_pdf
    - filt_cell_count_plot_png
    - filt_cell_count_plot_pdf
    - filt_umi_log10_density_plot_png
    - filt_umi_log10_density_plot_pdf
    - filt_gene_log10_density_plot_png
    - filt_gene_log10_density_plot_pdf
    - filt_gene_umi_log10_correlation_plot_png
    - filt_gene_umi_log10_correlation_plot_pdf
    - filt_mito_perc_log10_density_plot_png
    - filt_mito_perc_log10_density_plot_pdf
    - filt_novelty_score_log10_density_plot_png
    - filt_novelty_score_log10_density_plot_pdf
    - filt_qc_metrics_vln_plot_png
    - filt_qc_metrics_vln_plot_pdf
    - filt_qc_metrics_vln_plot_gr_by_cond_png
    - filt_qc_metrics_vln_plot_gr_by_cond_pdf
    - filt_unint_cell_cycle_eff_pca_plot_png
    - filt_unint_cell_cycle_eff_pca_plot_pdf
    - filt_unint_mito_perc_eff_pca_plot_png
    - filt_unint_mito_perc_eff_pca_plot_pdf
    - filt_unint_umap_plot_spl_by_ident_png
    - filt_unint_umap_plot_spl_by_ident_pdf
    - filt_int_elbow_plot_png
    - filt_int_elbow_plot_pdf
    - filt_int_pca_plot_png
    - filt_int_pca_plot_pdf
    - filt_int_pca_heatmap_png
    - filt_int_pca_heatmap_pdf
    - filt_int_pca_loadings_plot_png
    - filt_int_pca_loadings_plot_pdf
    - filt_int_umap_plot_spl_by_ident_png
    - filt_int_umap_plot_spl_by_ident_pdf
    - filt_int_cl_umap_plot_spl_by_cond_res_png
    - filt_int_cl_umap_plot_spl_by_cond_res_pdf
    - filt_int_cl_umap_ctype_pred_plot_res_png
    - filt_int_cl_umap_ctype_pred_plot_res_pdf
    - filt_int_cl_umap_plot_spl_by_ph_res_png
    - filt_int_cl_umap_plot_spl_by_ph_res_pdf
    - filt_int_cl_umap_qc_metrics_plot_res_png
    - filt_int_cl_umap_qc_metrics_plot_res_pdf
    - seurat_data_rds
    - conserved_gene_markers
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

s:name: "Seurat for comparative scRNA-seq analysis of across experimental conditions"
label: "Seurat for comparative scRNA-seq analysis of across experimental conditions"
s:alternateName: "Seurat for comparative scRNA-seq analysis of across experimental conditions"

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
  Runs Seurat for comparative scRNA-seq analysis of across experimental conditions
  ================================================================================