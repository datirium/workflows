cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/seurat:v0.0.1


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

  cell_cycle_data:
    type: File
    inputBinding:
      prefix: "--cellcycle"    
    doc: |
      Path to the TSV/CSV file with cell cycle data.
      First column - 'phase', second column 'gene_id'

  conditions_data:
    type: File
    inputBinding:
      prefix: "--condition"    
    doc: |
      Path to the TSV/CSV file to define datasets conditions
      for grouping. First column - 'library_id' with values
      from the --identity file, second column 'condition'.
      Default: each dataset is assigned to its own biological
      condition

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
      Include cells with the novelty score not lower that this
      value (calculated as log10(genes)/log10(UMIs)).
      Default: 0.8

  maximum_mito_perc:
    type: float?
    inputBinding:
      prefix: "--maxmt"
    doc: |
      Include cells with the mitochondrial contamination percentage
      not bigger that this value.
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

  only_positive_markers:
    type: boolean?
    inputBinding:
      prefix: "--onlypos"
    doc: |
      Return only positive markers when running conserved gene markers
      identification.
      Default: false

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
      glob: "*_raw_cell_count_plot.png"
    doc: |
      Raw number of cells per dataset plot in PNG format

  raw_cell_count_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_cell_count_plot.pdf"
    doc: |
      Raw number of cells per dataset plot in PDF format

  raw_umi_log10_density_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_umi_log10_density_plot.png"
    doc: |
      Raw UMI per cell log10 density plot in PNG format

  raw_umi_log10_density_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_umi_log10_density_plot.pdf"
    doc: |
      Raw UMI per cell log10 density plot in PDF format

  raw_gene_log10_density_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_log10_density_plot.png"
    doc: |
      Raw gene per cell log10 density plot in PNG format

  raw_gene_log10_density_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_log10_density_plot.pdf"
    doc: |
      Raw gene per cell log10 density plot in PDF format

  raw_gene_umi_log10_correlation_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi_log10_correlation_plot.png"
    doc: |
      Raw log10 gene vs log10 UMI correlation plot in PNG format

  raw_gene_umi_log10_correlation_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi_log10_correlation_plot.pdf"
    doc: |
      Raw log10 gene vs log10 UMI correlation plot in PDF format

  raw_mito_perc_log10_density_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_mito_perc_log10_density_plot.png"
    doc: |
      Raw mitochondrial gene percentage per cell log10 density plot in PNG format

  raw_mito_perc_log10_density_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_mito_perc_log10_density_plot.pdf"
    doc: |
      Raw mitochondrial gene percentage per cell log10 density plot in PDF format

  raw_novelty_score_log10_density_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_novelty_score_log10_density_plot.png"
    doc: |
      Raw novelty score (log10Gene/log10UMI) per cell log10 density plot in PNG format

  raw_novelty_score_log10_density_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_novelty_score_log10_density_plot.pdf"
    doc: |
      Raw novelty score (log10Gene/log10UMI) per cell log10 density plot in PDF format

  raw_qc_metrics_vln_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_qc_metrics_vln_plot.png"
    doc: |
      Raw QC metrics violin plot in PNG format

  raw_qc_metrics_vln_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_qc_metrics_vln_plot.pdf"
    doc: |
      Raw QC metrics violin plot in PDF format

  raw_qc_metrics_vln_plot_gr_by_cond_png:
    type: File?
    outputBinding:
      glob: "*_raw_qc_metrics_vln_plot_gr_by_cond.png"
    doc: |
      Raw QC metrics violin plot grouped by condition in PNG format

  raw_qc_metrics_vln_plot_gr_by_cond_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_qc_metrics_vln_plot_gr_by_cond.pdf"
    doc: |
      Raw QC metrics violin plot grouped by condition in PDF format


  filt_cell_count_plot_png:
    type: File?
    outputBinding:
      glob: "*_filt_cell_count_plot.png"
    doc: |
      Filtered number of cells per dataset plot in PNG format

  filt_cell_count_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_cell_count_plot.pdf"
    doc: |
      Filtered number of cells per dataset plot in PDF format

  filt_umi_log10_density_plot_png:
    type: File?
    outputBinding:
      glob: "*_filt_umi_log10_density_plot.png"
    doc: |
      Filtered UMI per cell log10 density plot in PNG format

  filt_umi_log10_density_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_umi_log10_density_plot.pdf"
    doc: |
      Filtered UMI per cell log10 density plot in PDF format

  filt_gene_log10_density_plot_png:
    type: File?
    outputBinding:
      glob: "*_filt_gene_log10_density_plot.png"
    doc: |
      Filtered gene per cell log10 density plot in PNG format

  filt_gene_log10_density_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_gene_log10_density_plot.pdf"
    doc: |
      Filtered gene per cell log10 density plot in PDF format

  filt_gene_umi_log10_correlation_plot_png:
    type: File?
    outputBinding:
      glob: "*_filt_gene_umi_log10_correlation_plot.png"
    doc: |
      Filtered log10 gene vs log10 UMI correlation plot in PNG format

  filt_gene_umi_log10_correlation_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_gene_umi_log10_correlation_plot.pdf"
    doc: |
      Filtered log10 gene vs log10 UMI correlation plot in PDF format

  filt_mito_perc_log10_density_plot_png:
    type: File?
    outputBinding:
      glob: "*_filt_mito_perc_log10_density_plot.png"
    doc: |
      Filtered mitochondrial gene percentage per cell log10 density plot in PNG format

  filt_mito_perc_log10_density_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_mito_perc_log10_density_plot.pdf"
    doc: |
      Filtered mitochondrial gene percentage per cell log10 density plot in PDF format

  filt_novelty_score_log10_density_plot_png:
    type: File?
    outputBinding:
      glob: "*_filt_novelty_score_log10_density_plot.png"
    doc: |
      Filtered novelty score (log10Gene/log10UMI) per cell log10 density plot in PNG format

  filt_novelty_score_log10_density_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_novelty_score_log10_density_plot.pdf"
    doc: |
      Filtered novelty score (log10Gene/log10UMI) per cell log10 density plot in PDF format

  filt_qc_metrics_vln_plot_png:
    type: File?
    outputBinding:
      glob: "*_filt_qc_metrics_vln_plot.png"
    doc: |
      Filtered QC metrics violin plot in PNG format

  filt_qc_metrics_vln_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_qc_metrics_vln_plot.pdf"
    doc: |
      Filtered QC metrics violin plot in PDF format

  filt_qc_metrics_vln_plot_gr_by_cond_png:
    type: File?
    outputBinding:
      glob: "*_filt_qc_metrics_vln_plot_gr_by_cond.png"
    doc: |
      Filtered QC metrics violin plot grouped by condition in PNG format

  filt_qc_metrics_vln_plot_gr_by_cond_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_qc_metrics_vln_plot_gr_by_cond.pdf"
    doc: |
      Filtered QC metrics violin plot grouped by condition in PDF format


  filt_unint_cell_cycle_eff_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_filt_unint_cell_cycle_eff_pca_plot.png"
    doc: |
      Filtered unintegrated PCA plot to evaluate cell cycle as a source of
      unwanted variation in PNG format

  filt_unint_cell_cycle_eff_pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_unint_cell_cycle_eff_pca_plot.pdf"
    doc: |
      Filtered unintegrated PCA plot to evaluate cell cycle as a source of
      unwanted variation in PDF format

  filt_unint_mito_perc_eff_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_filt_unint_mito_perc_eff_pca_plot.png"
    doc: |
      Filtered unintegrated PCA plot to evaluate mitochondrial contamination
      percentage as a source of unwanted variation in PNG format

  filt_unint_mito_perc_eff_pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_unint_mito_perc_eff_pca_plot.pdf"
    doc: |
      Filtered unintegrated PCA plot to evaluate mitochondrial contamination
      percentage as a source of unwanted variation in PDF format

  filt_unint_umap_plot_spl_by_ident_png:
    type: File?
    outputBinding:
      glob: "*_filt_unint_umap_plot_spl_by_ident.png"
    doc: |
      Filtered unintegrated UMAP plot split by cell identity in PNG format

  filt_unint_umap_plot_spl_by_ident_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_unint_umap_plot_spl_by_ident.pdf"
    doc: |
      Filtered unintegrated UMAP plot split by cell identity in PDF format


  filt_int_elbow_plot_png:
    type: File?
    outputBinding:
      glob: "*_filt_int_elbow_plot.png"
    doc: |
      Filtered integrated Elbow plot to evaluate data dimensionality
      in PNG format

  filt_int_elbow_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_int_elbow_plot.pdf"
    doc: |
      Filtered integrated Elbow plot to evaluate data dimensionality
      in PDF format

  filt_int_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_filt_int_pca_plot.png"    
    doc: |
      Filtered integrated PCA plot in PNG format

  filt_int_pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_int_pca_plot.pdf"    
    doc: |
      Filtered integrated PCA plot in PDF format

  filt_int_pca_heatmap_png:
    type: File?
    outputBinding:
      glob: "*_filt_int_pca_heatmap.png"        
    doc: |
      Filtered integrated PCA heatmap to evaluate data dimensionality
      in PNG format

  filt_int_pca_heatmap_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_int_pca_heatmap.pdf"        
    doc: |
      Filtered integrated PCA heatmap to evaluate data dimensionality
      in PDF format

  filt_int_pca_loadings_plot_png:
    type: File?
    outputBinding:
      glob: "*_filt_int_pca_loadings_plot.png"   
    doc: |
      Filtered integrated PCA loadings plot to evaluate data dimensionality
      in PNG format

  filt_int_pca_loadings_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_int_pca_loadings_plot.pdf"
    doc: |
      Filtered integrated PCA loadings plot to evaluate data dimensionality
      in PDF format

  filt_int_umap_plot_spl_by_ident_png:
    type: File?
    outputBinding:
      glob: "*_filt_int_umap_plot_spl_by_ident.png"
    doc: |
      Filtered integrated UMAP plot in PNG format

  filt_int_umap_plot_spl_by_ident_pdf:
    type: File?
    outputBinding:
      glob: "*_filt_int_umap_plot_spl_by_ident.pdf"
    doc: |
      Filtered integrated UMAP plot in PDF format


  filt_int_cl_umap_plot_spl_by_cond_res_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_filt_int_cl_umap_plot_spl_by_cond_res_*.png"
    doc: |
      Filtered integrated clustered UMAP plots with variable resolution
      in PNG format

  filt_int_cl_umap_plot_spl_by_cond_res_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_filt_int_cl_umap_plot_spl_by_cond_res_*.pdf"
    doc: |
      Filtered integrated clustered UMAP plots with variable resolution
      in PDF format

  filt_int_cl_umap_plot_spl_by_ph_res_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_filt_int_cl_umap_plot_spl_by_ph_res_*.png"
    doc: |
      Filtered integrated clustered UMAP plots split by cell cycle phase
      with variable resolution in PNG format

  filt_int_cl_umap_plot_spl_by_ph_res_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_filt_int_cl_umap_plot_spl_by_ph_res_*.pdf"
    doc: |
      Filtered integrated clustered UMAP plots split by cell cycle phase
      with variable resolution in PDF format

  filt_int_cl_umap_qc_metrics_plot_res_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_filt_int_cl_umap_qc_metrics_plot_res_*.png"
    doc: |
      Filtered integrated clustered UMAP QC metrics plots with variable
      resolution in PNG format

  filt_int_cl_umap_qc_metrics_plot_res_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_filt_int_cl_umap_qc_metrics_plot_res_*.pdf"
    doc: |
      Filtered integrated clustered UMAP QC metrics plots with variable
      resolution in PDF format

  
  seurat_data_rds:
    type: File?
    outputBinding:
      glob: "*_data.rds"
    doc: |
      Filtered integrated clustered Seurat data in RDS format

  conserved_gene_markers:
    type: File
    outputBinding:
      glob: "*_conserved_gene_markers.tsv"
    doc: |
      Conserved gene markers file for all clusters and all resolutions
      irrespective of condition in TSV format

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
      HTML index file from the directory with UCSC Cellbrowser
      formatted html data


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
  Runs Seurat for comparative scRNA-seq analysis of across experimental conditions
