cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/seurat-wnn:v0.0.1


inputs:

  feature_bc_matrices_folder:
    type: Directory
    inputBinding:
      prefix: "--mex"
    doc: |
      Path to the folder with feature-barcode matrix from Cell Ranger ARC Count
      in MEX format. The rows consist of all the gene and peak features concatenated
      together and the columns are restricted to those barcodes that are identified
      as cells.

  atac_fragments_file:
    type: File
    secondaryFiles:
    - .tbi
    inputBinding:
      prefix: "--fragments"
    doc: |
      Count and barcode information for every ATAC fragment observed in
      the experiment in TSV format. Tbi-index file is required.

  annotation_gtf_file:
    type: File
    inputBinding:
      prefix: "--annotations"
    doc: |
      Path to the genome annotation file in GTF format

  blacklisted_regions_file:
    type: File?
    inputBinding:
      prefix: "--blacklisted"
    doc: |
      Path to the blacklisted regions file in BED format

  gex_minimum_cells:
    type: int?
    inputBinding:
      prefix: "--gexmincells"
    doc: |
      Include only GEX features detected in at least this many cells.
      Default: 5

  gex_minimum_features:
    type: int?
    inputBinding:
      prefix: "--mingenes"
    doc: |
      Include cells where at least this many GEX features are detected.
      Default: 250

  gex_maximum_features:
    type: int?
    inputBinding:
      prefix: "--maxgenes"
    doc: |
      Include cells with the number of GEX features not bigger than this value.
      Default: 5000

  gex_minimum_umis:
    type: int?
    inputBinding:
      prefix: "--gexminumi"
    doc: |
      Include cells where at least this many GEX UMIs (transcripts) are detected.
      Default: 500

  mito_pattern:
    type: string?
    inputBinding:
      prefix: "--mitopattern"
    doc: |
      Regex pattern to identify mitochondrial GEX features.
      Default: '^Mt-'

  maximum_mito_perc:
    type: float?
    inputBinding:
      prefix: "--maxmt"
    doc: |
      Include cells with the percentage of GEX transcripts mapped to mitochondrial
      genes not bigger than this value.
      Default: 5

  minimum_novelty_score:
    type: float?
    inputBinding:
      prefix: "--minnovelty"
    doc: |
      Include cells with the novelty score not lower than this value
      calculated for GEX as log10(genes)/log10(UMIs).
      Default: 0.8

  atac_minimum_cells:
    type: int?
    inputBinding:
      prefix: "--atacmincells"
    doc: |
      Include only ATAC features detected in at least this many cells.
      Default: 5

  atac_minimum_umis:
    type: int?
    inputBinding:
      prefix: "--atacminumi"
    doc: |
      Include cells where at least this many ATAC UMIs (transcripts) are detected.
      Default: 1000

  maximum_nucl_signal:
    type: float?
    inputBinding:
      prefix: "--maxnuclsignal"
    doc: |
      Include cells with the nucleosome signal not bigger than this value.
      Nucleosome signal quantifies the approximate ratio of mononucleosomal
      to nucleosome-free fragments.
      Default: 4

  minimum_frip:
    type: float?
    inputBinding:
      prefix: "--minfrip"
    doc: |
      Include cells with the FRiP not lower than this value.
      Default: 0.15

  maximum_blacklisted_ratio:
    type: float?
    inputBinding:
      prefix: "--maxblacklisted"
    doc: |
      Include cells with the ratio of reads in genomic blacklist regions
      not bigger than this value.
      Default: 0.05

  call_peaks:
    type: boolean?
    inputBinding:
      prefix: "--callpeaks"
    doc: |
      Call peaks with MACS2 instead of those that are provided by Cell Ranger ARC Count.
      Default: false

  gex_selected_features:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--gexfeatures"
    doc: |
      GEX features of interest to evaluate expression.
      Default: None

  gex_high_var_features_count:
    type: int?
    inputBinding:
      prefix: "--highvarcount"
    doc: |
      Number of highly variable features to detect. Used for datasets integration,
      scaling, and dimensional reduction.
      Default: 3000

  gex_dimensionality:
    type: int?
    inputBinding:
      prefix: "--gexndim"
    doc: |
      Number of principal components to use in GEX UMAP projection and clustering
      (from 1 to 50). Use Elbow plot to adjust this parameter.
      Default: 50

  atac_dimensionality:
    type: int?
    inputBinding:
      prefix: "--atacndim"
    doc: |
      Number of principal components to use in ATAC UMAP projection and clustering
      (from 1 to 50). Use Elbow plot to adjust this parameter.
      Default: 50

  resolution:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--resolution"
    doc: |
      Clustering resolution. Can be set as an array.
      Default: 0.3

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

  raw_gex_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gex_umi_dnst.png"
    doc: |
      GEX UMI density per cell (not filtered).
      PNG format

  raw_gex_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gex_umi_dnst.pdf"
    doc: |
      GEX UMI density per cell (not filtered).
      PDF format

  raw_atac_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_atac_umi_dnst.png"
    doc: |
      ATAC UMI density per cell (not filtered).
      PNG format

  raw_atac_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_atac_umi_dnst.pdf"
    doc: |
      ATAC UMI density per cell (not filtered).
      PDF format

  raw_gene_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst.png"
    doc: |
      Gene density per cell (not filtered).
      PNG format

  raw_gene_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst.pdf"
    doc: |
      Gene density per cell (not filtered).
      PDF format
  
  raw_peak_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_peak_dnst.png"
    doc: |
      Peak density per cell (not filtered).
      PNG format

  raw_peak_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_peak_dnst.pdf"
    doc: |
      Peak density per cell (not filtered).
      PDF format

  raw_bl_cnts_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_bl_cnts_dnst.png"
    doc: |
      Density of fraction of reads within blacklisted regions per cell (not filtered).
      PNG format

  raw_bl_cnts_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_bl_cnts_dnst.pdf"
    doc: |
      Density of fraction of reads within blacklisted regions per cell (not filtered).
      PDF format

  raw_gex_atac_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gex_atac_umi_corr.png"
    doc: |
      GEX vs ATAC UMIs per cell correlation (not filtered).
      PNG format

  raw_gex_atac_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gex_atac_umi_corr.pdf"
    doc: |
      GEX vs ATAC UMIs per cell correlation (not filtered).
      PDF format

  raw_frg_len_hist_png:
    type: File?
    outputBinding:
      glob: "*_raw_frg_len_hist.png"
    doc: |
      Fragments Length Histogram (not filtered).
      PNG format

  raw_frg_len_hist_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_frg_len_hist.pdf"
    doc: |
      Fragments Length Histogram (not filtered).
      PDF format

  raw_gene_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi_corr.png"
    doc: |
      Genes vs GEX UMIs per cell correlation (not filtered).
      PNG format

  raw_gene_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi_corr.pdf"
    doc: |
      Genes vs GEX UMIs per cell correlation (not filtered).
      PDF format

  raw_mito_perc_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_mito_perc_dnst.png"
    doc: |
      Density of transcripts mapped to mitochondrial genes per cell (not filtered).
      PNG format

  raw_mito_perc_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_mito_perc_dnst.pdf"
    doc: |
      Density of transcripts mapped to mitochondrial genes per cell (not filtered).
      PDF format

  raw_nvlt_score_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_score_dnst.png"
    doc: |
      Novelty score density per cell (not filtered).
      PNG format

  raw_nvlt_score_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_score_dnst.pdf"
    doc: |
      Novelty score density per cell (not filtered).
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

  raw_qc_mtrcs_tsv:
    type: File?
    outputBinding:
      glob: "*_raw_qc_mtrcs.tsv"
    doc: |
      QC metrics densities per cell (not filtered).
      TSV format

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

  fltr_gex_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_gex_umi_dnst.png"
    doc: |
      GEX UMI density per cell (filtered).
      PNG format

  fltr_gex_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_gex_umi_dnst.pdf"
    doc: |
      GEX UMI density per cell (filtered).
      PDF format

  fltr_atac_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_atac_umi_dnst.png"
    doc: |
      ATAC UMI density per cell (filtered).
      PNG format

  fltr_atac_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_atac_umi_dnst.pdf"
    doc: |
      ATAC UMI density per cell (filtered).
      PDF format

  fltr_gene_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_dnst.png"
    doc: |
      Gene density per cell (filtered).
      PNG format

  fltr_gene_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_dnst.pdf"
    doc: |
      Gene density per cell (filtered).
      PDF format
  
  fltr_peak_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_peak_dnst.png"
    doc: |
      Peak density per cell (filtered).
      PNG format

  fltr_peak_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_peak_dnst.pdf"
    doc: |
      Peak density per cell (filtered).
      PDF format

  fltr_bl_cnts_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_bl_cnts_dnst.png"
    doc: |
      Density of fraction of reads within blacklisted regions per cell (filtered).
      PNG format

  fltr_bl_cnts_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_bl_cnts_dnst.pdf"
    doc: |
      Density of fraction of reads within blacklisted regions per cell (filtered).
      PDF format

  fltr_gex_atac_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_gex_atac_umi_corr.png"
    doc: |
      GEX vs ATAC UMIs per cell correlation (filtered).
      PNG format

  fltr_gex_atac_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_gex_atac_umi_corr.pdf"
    doc: |
      GEX vs ATAC UMIs per cell correlation (filtered).
      PDF format

  fltr_frg_len_hist_png:
    type: File?
    outputBinding:
      glob: "*_fltr_frg_len_hist.png"
    doc: |
      Fragments Length Histogram (filtered).
      PNG format

  fltr_frg_len_hist_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_frg_len_hist.pdf"
    doc: |
      Fragments Length Histogram (filtered).
      PDF format

  fltr_gene_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_umi_corr.png"
    doc: |
      Genes vs GEX UMIs per cell correlation (filtered).
      PNG format

  fltr_gene_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_umi_corr.pdf"
    doc: |
      Genes vs GEX UMIs per cell correlation (filtered).
      PDF format

  fltr_mito_perc_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_mito_perc_dnst.png"
    doc: |
      Density of transcripts mapped to mitochondrial genes per cell (filtered).
      PNG format

  fltr_mito_perc_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_mito_perc_dnst.pdf"
    doc: |
      Density of transcripts mapped to mitochondrial genes per cell (filtered).
      PDF format

  fltr_nvlt_score_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_nvlt_score_dnst.png"
    doc: |
      Novelty score density per cell (filtered).
      PNG format

  fltr_nvlt_score_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_nvlt_score_dnst.pdf"
    doc: |
      Novelty score density per cell (filtered).
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

  fltr_qc_mtrcs_tsv:
    type: File?
    outputBinding:
      glob: "*_fltr_qc_mtrcs.tsv"
    doc: |
      QC metrics densities per cell (filtered).
      TSV format

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

  clst_gex_umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_gex_umap_res_*.png"
    doc: |
      Clustered UMAP projected PCA of filtered GEX datasets.
      PNG format

  clst_gex_umap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_gex_umap_res_*.pdf"
    doc: |
      Clustered UMAP projected PCA of filtered GEX datasets.
      PDF format

  clst_atac_umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_atac_umap_res_*.png"
    doc: |
      Clustered UMAP projected LSI of filtered ATAC datasets.
      PNG format

  clst_atac_umap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_atac_umap_res_*.pdf"
    doc: |
      Clustered UMAP projected LSI of filtered ATAC datasets.
      PDF format

  clst_wnn_umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_wnn_umap_res_*.png"
    doc: |
      Clustered UMAP projected WNN.
      PNG format

  clst_wnn_umap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_clst_wnn_umap_res_*.pdf"
    doc: |
      Clustered UMAP projected WNN.
      PDF format

  seurat_clst_data_rds:
    type: File?
    outputBinding:
      glob: "*_clst_data.rds"
    doc: |
      Clustered filtered integrated/scaled Seurat data.
      RDS format

  expr_avg_per_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_avg_per_clst_res_*.png"
    doc: |
      Scaled average log normalized gene expression per cluster.
      PNG format

  expr_avg_per_clst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_avg_per_clst_res_*.pdf"
    doc: |
      Scaled average log normalized gene expression per cluster.
      PDF format

  expr_per_clst_cell_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_per_clst_cell_res_*.png"
    doc: |
      Log normalized gene expression per cell of clustered datasets.
      PNG format

  expr_per_clst_cell_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_per_clst_cell_res_*.pdf"
    doc: |
      Log normalized gene expression per cell of clustered datasets.
      PDF format

  expr_dnst_per_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_dnst_per_clst_res_*.png"
    doc: |
      Log normalized gene expression densities per cluster.
      PNG format

  expr_dnst_per_clst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_expr_dnst_per_clst_res_*.pdf"
    doc: |
      Log normalized gene expression densities per cluster.
      PDF format

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


baseCommand: ["run_seurat_wnn.R"]

stdout: seurat_wnn_stdout.log
stderr: seurat_wnn_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Seurat WNN Analysis"
s:name: "Seurat WNN Analysis"
s:alternateName: "Runs Seurat Weighted Nearest Neighbor Analysis"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/seurat-wnn-cluster.cwl
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
  Seurat WNN Analysis
  ===================
  Runs Seurat Weighted Nearest Neighbor Analysis


s:about: |
  usage: run_seurat_wnn.R [-h] --mex MEX --fragments FRAGMENTS
                                        --annotations ANNOTATIONS
                                        [--blacklisted BLACKLISTED]
                                        [--gexmincells GEXMINCELLS]
                                        [--mingenes MINGENES]
                                        [--maxgenes MAXGENES]
                                        [--gexminumi GEXMINUMI]
                                        [--mitopattern MITOPATTERN]
                                        [--maxmt MAXMT]
                                        [--minnovelty MINNOVELTY]
                                        [--atacmincells ATACMINCELLS]
                                        [--atacminumi ATACMINUMI]
                                        [--maxnuclsignal MAXNUCLSIGNAL]
                                        [--minfrip MINFRIP]
                                        [--maxblacklisted MAXBLACKLISTED]
                                        [--callpeaks]
                                        [--gexfeatures [GEXFEATURES [GEXFEATURES ...]]]
                                        [--highvarcount HIGHVARCOUNT]
                                        [--gexndim GEXNDIM]
                                        [--atacndim ATACNDIM]
                                        [--resolution [RESOLUTION [RESOLUTION ...]]]
                                        [--pdf] [--rds] [--output OUTPUT]
                                        [--threads THREADS]

  Runs Seurat Weighted Nearest Neighbor Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --mex MEX             Path to the folder with feature-barcode matrix from
                          Cell Ranger ARC Count in MEX format. The rows consist
                          of all the gene and peak features concatenated
                          together and the columns are restricted to those
                          barcodes that are identified as cells.
    --fragments FRAGMENTS
                          Count and barcode information for every ATAC fragment
                          observed in the experiment in TSV format. Tbi-index
                          file is required.
    --annotations ANNOTATIONS
                          Path to the genome annotation file in GTF format
    --blacklisted BLACKLISTED
                          Path to the blacklisted regions file in BED format
    --gexmincells GEXMINCELLS
                          Include only GEX features detected in at least this
                          many cells. Default: 5
    --mingenes MINGENES   Include cells where at least this many GEX features
                          are detected Default: 250
    --maxgenes MAXGENES   Include cells with the number of GEX features not
                          bigger than this value. Default: 5000
    --gexminumi GEXMINUMI
                          Include cells where at least this many GEX UMIs
                          (transcripts) are detected. Default: 500
    --mitopattern MITOPATTERN
                          Regex pattern to identify mitochondrial GEX features.
                          Default: '^Mt-'
    --maxmt MAXMT         Include cells with the percentage of GEX transcripts
                          mapped to mitochondrial genes not bigger than this
                          value. Default: 5
    --minnovelty MINNOVELTY
                          Include cells with the novelty score not lower than
                          this value, calculated for GEX as
                          log10(genes)/log10(UMIs). Default: 0.8
    --atacmincells ATACMINCELLS
                          Include only ATAC features detected in at least this
                          many cells. Default: 5
    --atacminumi ATACMINUMI
                          Include cells where at least this many ATAC UMIs
                          (transcripts) are detected. Default: 1000
    --maxnuclsignal MAXNUCLSIGNAL
                          Include cells with the nucleosome signal not bigger
                          than this value. Nucleosome signal quantifies the
                          approximate ratio of mononucleosomal to nucleosome-
                          free fragments. Default: 4
    --minfrip MINFRIP     Include cells with the FRiP not lower than this value.
                          Default: 0.15
    --maxblacklisted MAXBLACKLISTED
                          Include cells with the ratio of reads in genomic
                          blacklist regions not bigger than this value. Default:
                          0.05
    --callpeaks           Call peaks with MACS2 instead of those that are
                          provided by Cell Ranger ARC Count. Default: false
    --gexfeatures [GEXFEATURES [GEXFEATURES ...]]
                          GEX features of interest to evaluate expression.
                          Default: None
    --highvarcount HIGHVARCOUNT
                          Number of highly variable features to detect. Used for
                          datasets integration, scaling, and dimensional
                          reduction. Default: 3000
    --gexndim GEXNDIM     Number of principal components to use in GEX UMAP
                          projection and clustering (from 1 to 50). Use Elbow
                          plot to adjust this parameter. Default: 50
    --atacndim ATACNDIM   Number of principal components to use in ATAC UMAP
                          projection and clustering (from 1 to 50). Use Elbow
                          plot to adjust this parameter. Default: 50
    --resolution [RESOLUTION [RESOLUTION ...]]
                          Clustering resolution. Can be set as an array.
                          Default: 0.3
    --pdf                 Export plots in PDF. Default: false
    --rds                 Save Seurat data to RDS file. Default: false
    --output OUTPUT       Output prefix. Default: ./seurat
    --threads THREADS     Threads. Default: 1