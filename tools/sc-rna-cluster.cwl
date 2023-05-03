cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.20


inputs:

  query_data_rds:
    type: File
    inputBinding:
      prefix: "--query"
    doc: |
      Path to the RDS file to load Seurat object from. This file should include genes
      expression information stored in the RNA assay, as well as 'pca' and 'rnaumap'
      dimensionality reductions applied to that assay.

  dimensions:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--dimensions"
    doc: |
      Dimensionality to use when constructing nearest-
      neighbor graph before clustering (from 1 to 50). If
      single value N is provided, use from 1 to N
      dimensions. If multiple values are provided, subset to
      only selected dimensions.
      Default: from 1 to 10

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
      Distance metric used when constructing nearest-neighbor graph before clustering.
      Default: euclidean

  cluster_algorithm:
    type:
    - "null"
    - type: enum
      symbols:
      - "louvain"
      - "mult-louvain"
      - "slm"
      - "leiden"
    inputBinding:
      prefix: "--algorithm"
    doc: |
      Algorithm for modularity optimization when running clustering.
      Default: louvain

  resolution:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--resolution"
    doc: |
      Clustering resolution applied to the constructed nearest-neighbor graph.
      Can be set as an array but only the first item from the list will be used
      for cluster labels and gene markers in the UCSC Cell Browser when running
      with --cbbuild and --diffgenes parameters.
      Default: 0.3, 0.5, 1.0

  genes_of_interest:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--genes"
    doc: |
      Genes of interest to build genes expression plots.
      Default: None

  identify_diff_genes:
    type: boolean?
    inputBinding:
      prefix: "--diffgenes"
    doc: |
      Identify differentially expressed genes (putative gene markers) between each
      pair of clusters for all resolutions.
      Default: false

  minimum_logfc:
    type: float?
    inputBinding:
      prefix: "--logfc"
    doc: |
      For putative gene markers identification include only those genes that
      on average have log fold change difference in expression between every
      tested pair of clusters not lower than this value. Ignored if '--diffgenes'
      is not set.
      Default: 0.25

  minimum_pct:
    type: float?
    inputBinding:
      prefix: "--minpct"
    doc: |
      For putative gene markers identification include only those genes that
      are detected in not lower than this fraction of cells in either of the
      two tested clusters. Ignored if '--diffgenes' is not set.
      Default: 0.1

  only_positive_diff_genes:
    type: boolean?
    inputBinding:
      prefix: "--onlypos"
    doc: |
      For putative gene markers identification return only positive markers.
      Ignored if '--diffgenes' is not set.
      Default: false

  test_to_use:
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
      Statistical test to use for putative gene markers identification.
      Ignored if '--diffgenes' is not set.
      Default: wilcox

  export_pdf_plots:
    type: boolean?
    inputBinding:
      prefix: "--pdf"
    doc: |
      Export plots in PDF.
      Default: false

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
    inputBinding:
      prefix: "--theme"
    doc: |
      Color theme for all generated plots. One of gray, bw, linedraw, light,
      dark, minimal, classic, void.
      Default: classic

  verbose:
    type: boolean?
    inputBinding:
      prefix: "--verbose"
    doc: |
      Print debug information.
      Default: false

  export_h5seurat_data:
    type: boolean?
    inputBinding:
      prefix: "--h5seurat"
    doc: |
      Save Seurat data to h5seurat file.
      Default: false

  export_h5ad_data:
    type: boolean?
    inputBinding:
      prefix: "--h5ad"
    doc: |
      Save Seurat data to h5ad file.
      Default: false

  export_scope_data:
    type: boolean?
    inputBinding:
      prefix: "--scope"
    doc: |
      Save Seurat data to SCope compatible loom file.
      Default: false

  export_ucsc_cb:
    type: boolean?
    inputBinding:
      prefix: "--cbbuild"
    doc: |
      Export results to UCSC Cell Browser. Default: false

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output prefix.
      Default: ./sc

  parallel_memory_limit:
    type: int?
    inputBinding:
      prefix: "--memory"
    doc: |
      Maximum memory in GB allowed to be shared between the workers
      when using multiple --cpus.
      Default: 32

  vector_memory_limit:
    type: int?
    default: 128
    doc: |
      Maximum vector memory in GB allowed to be used by R.
      Default: 128

  threads:
    type: int?
    inputBinding:
      prefix: "--cpus"
    doc: |
      Number of cores/cpus to use.
      Default: 1


outputs:

  umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_res_*.png"
    doc: |
      Clustered cells UMAP.
      PNG format

  umap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_res_*.pdf"
    doc: |
      Clustered cells UMAP.
      PDF format

  slh_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_slh_res_*.png"
    doc: |
      Silhouette scores. Downsampled to max 500 cells per cluster.
      PNG format

  slh_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_slh_res_*.pdf"
    doc: |
      Silhouette scores. Downsampled to max 500 cells per cluster.
      PDF format

  umap_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_idnt_res_*.png"
    doc: |
      Split by dataset clustered cells UMAP.
      PNG format

  umap_spl_idnt_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_idnt_res_*.pdf"
    doc: |
      Split by dataset clustered cells UMAP.
      PDF format

  cmp_gr_clst_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_clst_spl_idnt_res_*.png"
    doc: |
      Grouped by cluster split by dataset cells composition plot. Downsampled.
      PNG format

  cmp_gr_clst_spl_idnt_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_clst_spl_idnt_res_*.pdf"
    doc: |
      Grouped by cluster split by dataset cells composition plot. Downsampled.
      PDF format

  cmp_gr_idnt_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_idnt_spl_clst_res_*.png"
    doc: |
      Grouped by dataset split by cluster cells composition plot. Downsampled.
      PNG format

  cmp_gr_idnt_spl_clst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_idnt_spl_clst_res_*.pdf"
    doc: |
      Grouped by dataset split by cluster cells composition plot. Downsampled.
      PDF format

  umap_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_cnd_res_*.png"
    doc: |
      Split by grouping condition clustered cells UMAP.
      PNG format

  umap_spl_cnd_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_cnd_res_*.pdf"
    doc: |
      Split by grouping condition clustered cells UMAP.
      PDF format

  cmp_gr_clst_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_clst_spl_cnd_res_*.png"
    doc: |
      Grouped by cluster split by condition cells composition plot. Downsampled.
      PNG format

  cmp_gr_clst_spl_cnd_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_clst_spl_cnd_res_*.pdf"
    doc: |
      Grouped by cluster split by condition cells composition plot. Downsampled.
      PDF format

  cmp_gr_cnd_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_cnd_spl_clst_res_*.png"
    doc: |
      Grouped by condition split by cluster cells composition plot. Downsampled.
      PNG format

  cmp_gr_cnd_spl_clst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_cnd_spl_clst_res_*.pdf"
    doc: |
      Grouped by condition split by cluster cells composition plot. Downsampled.
      PDF format

  umap_spl_ph_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_ph_res_*.png"
    doc: |
      Split by cell cycle phase clustered cells UMAP.
      PNG format

  umap_spl_ph_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_ph_res_*.pdf"
    doc: |
      Split by cell cycle phase clustered cells UMAP.
      PDF format

  cmp_gr_ph_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_idnt.png"
    doc: |
      Grouped by cell cycle phase split by dataset cells composition plot. Downsampled.
      PNG format

  cmp_gr_ph_spl_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_idnt.pdf"
    doc: |
      Grouped by cell cycle phase split by dataset cells composition plot. Downsampled.
      PDF format

  cmp_gr_ph_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_ph_spl_clst_res_*.png"
    doc: |
      Grouped by cell cycle phase split by cluster cells composition plot. Downsampled.
      PNG format

  cmp_gr_ph_spl_clst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_ph_spl_clst_res_*.pdf"
    doc: |
      Grouped by cell cycle phase split by cluster cells composition plot. Downsampled.
      PDF format

  xpr_avg_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_avg_res_*.png"
    doc: |
      Log normalized scaled average gene expression per cluster.
      PNG format

  xpr_avg_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_avg_res_*.pdf"
    doc: |
      Log normalized scaled average gene expression per cluster.
      PDF format

  xpr_per_cell_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_[!sgnl_]*.png"
    doc: |
      Log normalized gene expression on cells UMAP.
      PNG format

  xpr_per_cell_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_[!sgnl_]*.pdf"
    doc: |
      Log normalized gene expression on cells UMAP.
      PDF format

  xpr_per_cell_sgnl_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_sgnl_*.png"
    doc: |
      Log normalized gene expression density on cells UMAP.
      PNG format

  xpr_per_cell_sgnl_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_sgnl_*.pdf"
    doc: |
      Log normalized gene expression density on cells UMAP.
      PDF format

  xpr_dnst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_dnst_res_*.png"
    doc: |
      Log normalized gene expression density per cluster.
      PNG format

  xpr_dnst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_dnst_res_*.pdf"
    doc: |
      Log normalized gene expression density per cluster.
      PDF format

  xpr_htmp_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_htmp_res_*.png"
    doc: |
      Normalized gene expression heatmap grouped by cluster.
      PNG format

  xpr_htmp_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_htmp_res_*.pdf"
    doc: |
      Normalized gene expression heatmap grouped by cluster.
      PDF format

  gene_markers_tsv:
    type: File?
    outputBinding:
      glob: "*_gene_markers.tsv"
    doc: |
      Differentially expressed genes between each pair of clusters for all resolutions.
      TSV format

  ucsc_cb_config_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser"
    doc: |
      Directory with UCSC Cellbrowser configuration data.

  ucsc_cb_html_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser/html_data"
    doc: |
      Directory with UCSC Cellbrowser html data.

  ucsc_cb_html_file:
    type: File?
    outputBinding:
      glob: "*_cellbrowser/html_data/index.html"
    doc: |
      HTML index file from the directory with UCSC Cellbrowser html data.

  seurat_data_rds:
    type: File
    outputBinding:
      glob: "*_data.rds"
    doc: |
      Reduced Seurat data in RDS format

  seurat_data_h5seurat:
    type: File?
    outputBinding:
      glob: "*_data.h5seurat"
    doc: |
      Reduced Seurat data in h5seurat format

  seurat_data_h5ad:
    type: File?
    outputBinding:
      glob: "*_data.h5ad"
    doc: |
      Reduced Seurat data in h5ad format

  seurat_data_scope:
    type: File?
    outputBinding:
      glob: "*_data.loom"
    doc: |
      Reduced Seurat data in SCope compatible loom format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sc_rna_cluster.R"]

stdout: sc_rna_cluster_stdout.log
stderr: sc_rna_cluster_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell RNA-Seq Cluster Analysis"
s:name: "Single-cell RNA-Seq Cluster Analysis"
s:alternateName: "Clusters single-cell RNA-Seq datasets, identifies gene markers"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-rna-cluster.cwl
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
  Single-cell RNA-Seq Cluster Analysis

  Clusters single-cell RNA-Seq datasets, identifies gene markers.


s:about: |
  usage: sc_rna_cluster.R
        [-h] --query QUERY [--dimensions [DIMENSIONS [DIMENSIONS ...]]]
        [--ametric {euclidean,cosine,manhattan,hamming}]
        [--algorithm {louvain,mult-louvain,slm,leiden}]
        [--resolution [RESOLUTION [RESOLUTION ...]]]
        [--genes [GENES [GENES ...]]] [--diffgenes] [--logfc LOGFC]
        [--minpct MINPCT] [--onlypos]
        [--testuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}]
        [--pdf] [--verbose] [--h5seurat] [--h5ad] [--cbbuild] [--output OUTPUT]
        [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
        [--cpus CPUS] [--memory MEMORY]

  Single-cell RNA-Seq Cluster Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include genes expression information
                          stored in the RNA assay, as well as 'pca' and
                          'rnaumap' dimensionality reductions applied to that
                          assay.
    --dimensions [DIMENSIONS [DIMENSIONS ...]]
                          Dimensionality to use when constructing nearest-
                          neighbor graph before clustering (from 1 to 50). If
                          single value N is provided, use from 1 to N
                          dimensions. If multiple values are provided, subset to
                          only selected dimensions. Default: from 1 to 10
    --ametric {euclidean,cosine,manhattan,hamming}
                          Distance metric used when constructing nearest-
                          neighbor graph before clustering. Default: euclidean
    --algorithm {louvain,mult-louvain,slm,leiden}
                          Algorithm for modularity optimization when running
                          clustering. Default: louvain
    --resolution [RESOLUTION [RESOLUTION ...]]
                          Clustering resolution applied to the constructed
                          nearest-neighbor graph. Can be set as an array but
                          only the first item from the list will be used for
                          cluster labels and gene markers in the UCSC Cell
                          Browser when running with --cbbuild and --diffgenes
                          parameters. Default: 0.3, 0.5, 1.0
    --genes [GENES [GENES ...]]
                          Genes of interest to build genes expression plots.
                          Default: None
    --diffgenes           Identify differentially expressed genes (putative gene
                          markers) between each pair of clusters for all
                          resolutions. Default: false
    --logfc LOGFC         For putative gene markers identification include only
                          those genes that on average have log fold change
                          difference in expression between every tested pair of
                          clusters not lower than this value. Ignored if '--
                          diffgenes' is not set. Default: 0.25
    --minpct MINPCT       For putative gene markers identification include only
                          those genes that are detected in not lower than this
                          fraction of cells in either of the two tested
                          clusters. Ignored if '--diffgenes' is not set.
                          Default: 0.1
    --onlypos             For putative gene markers identification return only
                          positive markers. Ignored if '--diffgenes' is not set.
                          Default: false
    --testuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}
                          Statistical test to use for putative gene markers
                          identification. Ignored if '--diffgenes' is not set.
                          Default: wilcox
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --h5seurat            Save Seurat data to h5seurat file. Default: false
    --h5ad                Save Seurat data to h5ad file. Default: false
    --cbbuild             Export results to UCSC Cell Browser. Default: false
    --output OUTPUT       Output prefix. Default: ./sc
    --theme {gray,bw,linedraw,light,dark,minimal,classic,void}
                          Color theme for all generated plots. Default: classic
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32