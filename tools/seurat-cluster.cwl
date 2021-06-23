cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/seurat:v0.0.4


inputs:

  feature_bc_matrices_folder:
    type: Directory
    inputBinding:
      prefix: "--mex"
    doc: |
      Path to the folder with not normalized aggregated
      feature-barcode matrices in MEX format

  aggregation_metadata:
    type: File
    inputBinding:
      prefix: "--identity"    
    doc: |
      Path to the aggregation CSV file to set the initial
      cell identity classes

  conditions_data:
    type: File?
    inputBinding:
      prefix: "--condition"    
    doc: |
      Path to the TSV/CSV file to define datasets conditions
      for grouping. First column - 'library_id' with values
      from the --identity file, second column 'condition'.
      Default: each dataset is assigned to its own biological
      condition

  classifier_rds:
    type: File?
    inputBinding:
      prefix: "--classifier"
    doc: |
      Path to the Garnett classifier rds file for cell type prediction.
      Default: skip cell type prediction

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
      Default: the same as "none" - do not convert gene names

  cell_cycle_data:
    type: File?
    inputBinding:
      prefix: "--cellcycle"
    doc: |
      Path to the TSV/CSV file with cell cycle data.
      First column - 'phase', second column 'gene_id'.
      Default: skip cell cycle score assignment

  barcodes_data:
    type: File?
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the headerless TSV/CSV file with selected barcodes
      (one per line) to prefilter input feature-barcode matrices.
      Default: use all cells

  minimum_cells:
    type: int?
    inputBinding:
      prefix: "--mincells"
    doc: |
      Include features detected in at least this many cells
      (applied to thoughout all datasets together).
      Default: 10

  minimum_features:
    type: int?
    inputBinding:
      prefix: "--minfeatures"
    doc: |
      Include cells where at least this many features are detected.
      Default: 250

  maximum_features:
    type: int?
    inputBinding:
      prefix: "--maxfeatures"
    doc: |
      Include cells with the number of features not bigger than this value.
      Default: 5000

  minimum_umis:
    type: int?
    inputBinding:
      prefix: "--minumi"
    doc: |
      Include cells where at least this many UMI are detected.
      Default: 500  

  minimum_novelty_score:
    type: float?
    inputBinding:
      prefix: "--minnovelty"
    doc: |
      Include cells with the novelty score not lower than this
      value (calculated as log10(genes)/log10(UMIs)).
      Default: 0.8

  maximum_mito_perc:
    type: float?
    inputBinding:
      prefix: "--maxmt"
    doc: |
      Include cells with the mitochondrial contamination percentage
      not bigger than this value.
      Default: 5

  mito_pattern:
    type: string?
    inputBinding:
      prefix: "--mitopattern"
    doc: |
      Regex pattern to identify mitochondrial reads.
      Default: ^Mt-

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
      Regress mitochondrial gene expression as a confounding source
      of variation.
      Default: false

  high_var_features_count:
    type: int?
    inputBinding:
      prefix: "--highvarcount"
    doc: |
      Number of higly variable features to detect.
      Default: 3000

  dimensionality:
    type: int?
    inputBinding:
      prefix: "--ndim"
    doc: |
      Number of principal components to use in clustering (1:50).
      Use Elbow plot to adjust this parameter.
      Default: 10

  resolution:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--resolution"
    doc: |
      Clustering resolution. Can be set as array.
      Default: 0.4 0.6 0.8 1.0 1.4

  minimum_logfc:
    type: float?
    inputBinding:
      prefix: "--logfc"
    doc: |
      Log fold change threshold for conserved gene markers identification.
      Default: 0.25

  minimum_pct:
    type: float?
    inputBinding:
      prefix: "--minpct"
    doc: |
      Minimum fraction of cells where genes used for conserved gene markers
      identification should be detected in either of two tested clusters.
      Default: 0.1

  only_positive_markers:
    type: boolean?
    inputBinding:
      prefix: "--onlypos"
    doc: |
      Return only positive markers when running conserved gene markers
      identification.
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
      Set test type to use for putative and conserved gene marker identification.
      Default: wilcox

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
      Split by condition mitochondrial gene percentage density per cell (not filtered).
      PNG format

  raw_mito_perc_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_mito_perc_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition mitochondrial gene percentage density per cell (not filtered).
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
      Split by condition mitochondrial gene percentage density per cell (filtered).
      PNG format

  fltr_mito_perc_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_mito_perc_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition mitochondrial gene percentage density per cell (filtered).
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
      Split by cell cycle phase PCA of filtered unintegrated datasets.
      PNG format

  fltr_pca_spl_by_ph_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_pca_spl_by_ph.pdf"
    doc: |
      Split by cell cycle phase PCA of filtered unintegrated datasets.
      PDF format

  fltr_pca_spl_by_mito_perc_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_pca_spl_by_mito_perc.png"
    doc: |
      Split by level of mitochondrial gene expression PCA of filtered unintegrated datasets.
      PNG format

  fltr_pca_spl_by_mito_perc_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_pca_spl_by_mito_perc.pdf"
    doc: |
      Split by level of mitochondrial gene expression PCA of filtered unintegrated datasets.
      PDF format

  fltr_umap_spl_by_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_umap_spl_by_idnt.png"
    doc: |
      Split by identity UMAP projected PCA of filtered unintegrated datasets.
      PNG format

  fltr_umap_spl_by_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_umap_spl_by_idnt.pdf"
    doc: |
      Split by identity UMAP projected PCA of filtered unintegrated datasets.
      PDF format


  ntgr_elbow_plot_png:
    type: File?
    outputBinding:
      glob: "*_ntgr_elbow.png"
    doc: |
      Elbow plot from PCA of filtered integrated datasets.
      PNG format

  ntgr_elbow_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_ntgr_elbow.pdf"
    doc: |
      Elbow plot from PCA of filtered integrated datasets.
      PDF format

  ntgr_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_ntgr_pca.png"    
    doc: |
      PCA of filtered integrated datasets.
      PNG format

  ntgr_pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_ntgr_pca.pdf"    
    doc: |
      PCA of filtered integrated datasets.
      PDF format

  ntgr_pca_heatmap_png:
    type: File?
    outputBinding:
      glob: "*_ntgr_pca_heatmap.png"        
    doc: |
      Genes per cells expression heatmap sorted by their PC scores from PCA of filtered integrated datasets.
      PNG format

  ntgr_pca_heatmap_pdf:
    type: File?
    outputBinding:
      glob: "*_ntgr_pca_heatmap.pdf"        
    doc: |
      Genes per cells expression heatmap sorted by their PC scores from PCA of filtered integrated datasets.
      PDF format

  ntgr_pca_loadings_plot_png:
    type: File?
    outputBinding:
      glob: "*_ntgr_pca_loadings.png"   
    doc: |
      PC scores of the most variant genes from PCA of filtered integrated datasets.
      PNG format

  ntgr_pca_loadings_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_ntgr_pca_loadings.pdf"
    doc: |
      PC scores of the most variant genes from PCA of filtered integrated datasets.
      PDF format

  ntgr_umap_spl_by_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_ntgr_umap_spl_by_idnt.png"
    doc: |
      Split by identity UMAP projected PCA of filtered integrated datasets.
      PNG format

  ntgr_umap_spl_by_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_ntgr_umap_spl_by_idnt.pdf"
    doc: |
      Split by identity UMAP projected PCA of filtered integrated datasets.
      PDF format


  clst_umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_res_*.png"
    doc: |
      Clustered UMAP projected PCA of filtered integrated datasets.
      PNG format

  clst_umap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_res_*.pdf"
    doc: |
      Clustered UMAP projected PCA of filtered integrated datasets.
      PDF format

  clst_umap_spl_by_cond_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_spl_by_cond_res_*.png"
    doc: |
      Split by condition clustered UMAP projected PCA of filtered integrated datasets.
      PNG format

  clst_umap_spl_by_cond_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_spl_by_cond_res_*.pdf"
    doc: |
      Split by condition clustered UMAP projected PCA of filtered integrated datasets.
      PDF format

  clst_umap_ctype_pred_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_ctype_pred_res_*.png"
    doc: |
      Grouped by predicted cell types UMAP projected PCA of filtered integrated datasets.
      PNG format

  clst_umap_ctype_pred_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_ctype_pred_res_*.pdf"
    doc: |
      Grouped by predicted cell types UMAP projected PCA of filtered integrated datasets.
      PDF format

  clst_umap_spl_by_ph_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_spl_by_ph_res_*.png"
    doc: |
      Split by cell cycle phase clustered UMAP projected PCA of filtered integrated datasets.
      PNG format

  clst_umap_spl_by_ph_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_umap_spl_by_ph_res_*.pdf"
    doc: |
      Split by cell cycle phase clustered UMAP projected PCA of filtered integrated datasets.
      PDF format

  clst_qc_mtrcs_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_qc_mtrcs_res_*.png"
    doc: |
      QC metrics for clustered UMAP projected PCA of filtered integrated datasets.
      PNG format

  clst_qc_mtrcs_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_qc_mtrcs_res_*.pdf"
    doc: |
      QC metrics for clustered UMAP projected PCA of filtered integrated datasets.
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
      Clustered filtered integrated Seurat data.
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


s:name: "seurat-cluster"
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
  Runs Seurat for comparative scRNA-seq analysis of across experimental conditions

 
s:about: |
  usage: run_seurat.R
        [-h] --mex MEX --identity IDENTITY [--condition CONDITION]
        [--classifier CLASSIFIER] [--cellcycle CELLCYCLE] [--barcodes BARCODES]
        [--mincells MINCELLS] [--minfeatures MINFEATURES]
        [--maxfeatures MAXFEATURES] [--minumi MINUMI] [--minnovelty MINNOVELTY]
        [--maxmt MAXMT] [--mitopattern MITOPATTERN] [--regresscellcycle]
        [--regressmt] [--highvarcount HIGHVARCOUNT] [--ndim NDIM]
        [--resolution [RESOLUTION [RESOLUTION ...]]] [--logfc LOGFC]
        [--minpct MINPCT] [--onlypos]
        [--testuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}]
        [--species {hs,mm,none}] [--pdf] [--rds] [--output OUTPUT]
        [--threads THREADS]

  Runs Seurat for comparative scRNA-seq analysis of across experimental
  conditions

  optional arguments:
    -h, --help            show this help message and exit
    --mex MEX             Path to the folder with not normalized aggregated
                          feature-barcode matrices in MEX format
    --identity IDENTITY   Path to the aggregation CSV file to set the initial
                          cell identity classes
    --condition CONDITION
                          Path to the TSV/CSV file to define datasets conditions
                          for grouping. First column - 'library_id' with values
                          from the --identity file, second column 'condition'.
                          Default: each dataset is assigned to its own
                          biological condition
    --classifier CLASSIFIER
                          Path to the Garnett classifier rds file for cell type
                          prediction. Default: skip cell type prediction
    --cellcycle CELLCYCLE
                          Path to the TSV/CSV file with cell cycle data. First
                          column - 'phase', second column 'gene_id'. Default:
                          skip cell cycle score assignment
    --barcodes BARCODES   Path to the headerless TSV/CSV file with selected
                          barcodes (one per line) to prefilter input feature-
                          barcode matrices. Default: use all cells
    --mincells MINCELLS   Include features detected in at least this many cells
                          (applied to thoughout all datasets together). Default:
                          10
    --minfeatures MINFEATURES
                          Include cells where at least this many features are
                          detected. Default: 250
    --maxfeatures MAXFEATURES
                          Include cells with the number of features not bigger
                          than this value. Default: 5000
    --minumi MINUMI       Include cells where at least this many UMI are
                          detected. Default: 500
    --minnovelty MINNOVELTY
                          Include cells with the novelty score not lower than
                          this value (calculated as log10(genes)/log10(UMIs)).
                          Default: 0.8
    --maxmt MAXMT         Include cells with the mitochondrial contamination
                          percentage not bigger than this value. Default: 5
    --mitopattern MITOPATTERN
                          Regex pattern to identify mitochondrial reads.
                          Default: ^Mt-
    --regresscellcycle    Regress cell cycle as a confounding source of
                          variation. Default: false
    --regressmt           Regress mitochondrial gene expression as a confounding
                          source of variation. Default: false
    --highvarcount HIGHVARCOUNT
                          Number of higly variable features to detect. Default:
                          3000
    --ndim NDIM           Number of principal components to use in clustering
                          (1:50). Use Elbow plot to adjust this parameter.
                          Default: 10
    --resolution [RESOLUTION [RESOLUTION ...]]
                          Clustering resolution. Can be set as array. Default:
                          0.4 0.6 0.8 1.0 1.4
    --logfc LOGFC         Log fold change threshold for conserved gene markers
                          identification. Default: 0.25
    --minpct MINPCT       Minimum fraction of cells where genes used for
                          conserved gene markers identification should be
                          detected in either of two tested clusters. Default:
                          0.1
    --onlypos             Return only positive markers when running conserved
                          gene markers identification. Default: false
    --testuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}
                          Set test type to use for putative and conserved gene
                          marker identification. Default: wilcox
    --species {hs,mm,none}
                          Select species for gene name conversion when running
                          cell type prediction with Garnett classifier. Default:
                          do not convert gene names
    --pdf                 Export plots in PDF. Default: false
    --rds                 Save Seurat data to RDS file. Default: false
    --output OUTPUT       Output prefix. Default: ./seurat
    --threads THREADS     Threads. Default: 1