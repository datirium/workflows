cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.29


inputs:

  query_data_rds:
    type: File
    inputBinding:
      prefix: "--query"
    doc: |
      Path to the RDS file to load Seurat object from. This file should include
      genes expression and chromatin accessibility information stored in the RNA
      and ATAC assays correspondingly. Additionally, 'pca', 'rnaumap', 'atac_lsi'
      and 'atacumap' dimensionality reductions should be present.

  rna_dimensions:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--rnadimensions"
    doc: |
      Dimensionality from the 'pca' reduction to use when constructing weighted
      nearest-neighbor graph before clustering (from 1 to 50). If single value N
      is provided, use from 1 to N dimensions. If multiple values are provided,
      subset to only selected dimensions.
      Default: from 1 to 10

  atac_dimensions:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--atacdimensions"
    doc: |
      Dimensionality from the 'atac_lsi' reduction to use when constructing weighted
      nearest-neighbor graph before clustering (from 1 to 50). If single value N
      is provided, use from 2 to N dimensions. If multiple values are provided,
      subset to only selected dimensions.
      Default: from 2 to 10

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
      Default: slm

  umap_spread:
    type: float?
    inputBinding:
      prefix: "--uspread"
    doc: |
      The effective scale of embedded points on UMAP. In combination with '--mindist'
      it determines how clustered/clumped the embedded points are.
      Default: 1

  umap_mindist:
    type: float?
    inputBinding:
      prefix: "--umindist"
    doc: |
      Controls how tightly the embedding is allowed compress points together on UMAP.
      Larger values ensure embedded points are moreevenly distributed, while smaller
      values allow the algorithm to optimise more accurately with regard to local structure.
      Sensible values are in the range 0.001 to 0.5.
      Default:  0.3

  umap_neighbors:
    type: int?
    inputBinding:
      prefix: "--uneighbors"
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
      UMAP implementation to run. If set to 'umap-learn' use --umetric 'correlation'
      Default: uwot

  resolution:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--resolution"
    doc: |
      Clustering resolution applied to the constructed weighted nearest-neighbor
      graph. Can be set as an array but only the first item from the list will
      be used for cluster labels and gene/peak markers in the UCSC Cell Browser
      when running with --cbbuild and --diffgenes/--diffpeaks parameters.
      Default: 0.3, 0.5, 1.0

  atac_fragments_file:
    type: File?
    secondaryFiles:
    - .tbi
    inputBinding:
      prefix: "--fragments"
    doc: |
      Count and barcode information for every ATAC fragment used in the loaded Seurat
      object. File should be saved in TSV format with tbi-index file.

  genes_of_interest:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--genes"
    doc: |
      Genes of interest to build gene expression and Tn5 insertion frequency plots
      for the nearest peaks. If '--fragments' is not provided only gene expression
      plots will be built.
      Default: None

  identify_diff_genes:
    type: boolean?
    inputBinding:
      prefix: "--diffgenes"
    doc: |
      Identify differentially expressed genes (putative gene markers) between each
      pair of clusters for all resolutions.
      Default: false

  identify_diff_peaks:
    type: boolean?
    inputBinding:
      prefix: "--diffpeaks"
    doc: |
      Identify differentially accessible peaks between each pair of clusters for all resolutions.
      Default: false

  rna_minimum_logfc:
    type: float?
    inputBinding:
      prefix: "--rnalogfc"
    doc: |
      For putative gene markers identification include only those genes that
      on average have log fold change difference in expression between every
      tested pair of clusters not lower than this value. Ignored if '--diffgenes'
      is not set.
      Default: 0.25

  rna_minimum_pct:
    type: float?
    inputBinding:
      prefix: "--rnaminpct"
    doc: |
      For putative gene markers identification include only those genes that
      are detected in not lower than this fraction of cells in either of the
      two tested clusters. Ignored if '--diffgenes' is not set.
      Default: 0.1

  only_positive_diff_genes:
    type: boolean?
    inputBinding:
      prefix: "--rnaonlypos"
    doc: |
      For putative gene markers identification return only positive markers.
      Ignored if '--diffgenes' is not set.
      Default: false

  rna_test_to_use:
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
      prefix: "--rnatestuse"
    doc: |
      Statistical test to use for putative gene markers identification.
      Ignored if '--diffgenes' is not set.
      Default: wilcox

  atac_minimum_logfc:
    type: float?
    inputBinding:
      prefix: "--ataclogfc"
    doc: |
      For differentially accessible peaks identification include only those peaks that
      on average have log fold change difference in the chromatin accessibility between
      every tested pair of clusters not lower than this value. Ignored if '--diffpeaks'
      is not set.
      Default: 0.25

  atac_minimum_pct:
    type: float?
    inputBinding:
      prefix: "--atacminpct"
    doc: |
      For differentially accessible peaks identification include only those peaks that
      are detected in not lower than this fraction of cells in either of the two tested
      clusters. Ignored if '--diffpeaks' is not set.
      Default: 0.05

  atac_test_to_use:
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
      prefix: "--atactestuse"
    doc: |
      Statistical test to use for differentially accessible peaks identification.
      Ignored if '--diffpeaks' is not set.
      Default: LR

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
      Only not normalized raw counts from the RNA assay
      will be saved. Default: false

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
      UMAP, colored by cluster.
      PNG format

  umap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_res_*.pdf"
    doc: |
      UMAP, colored by cluster.
      PDF format

  umap_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_idnt_res_*.png"
    doc: |
      UMAP, colored by cluster,
      split by dataset.
      PNG format

  umap_spl_idnt_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_idnt_res_*.pdf"
    doc: |
      UMAP, colored by cluster,
      split by dataset.
      PDF format

  cmp_gr_clst_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_clst_spl_idnt_res_*.png"
    doc: |
      Composition plot, colored by
      cluster, split by dataset,
      downsampled.
      PNG format

  cmp_gr_clst_spl_idnt_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_clst_spl_idnt_res_*.pdf"
    doc: |
      Composition plot, colored by
      cluster, split by dataset,
      downsampled.
      PDF format

  cmp_gr_idnt_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_idnt_spl_clst_res_*.png"
    doc: |
      Composition plot, colored by
      dataset, split by cluster,
      downsampled.
      PNG format

  cmp_gr_idnt_spl_clst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_idnt_spl_clst_res_*.pdf"
    doc: |
      Composition plot, colored by
      dataset, split by cluster,
      downsampled.
      PDF format

  umap_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_cnd_res_*.png"
    doc: |
      UMAP, colored by cluster,
      split by grouping condition.
      PNG format

  umap_spl_cnd_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_cnd_res_*.pdf"
    doc: |
      UMAP, colored by cluster,
      split by grouping condition.
      PDF format

  cmp_gr_clst_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_clst_spl_cnd_res_*.png"
    doc: |
      Composition plot, colored by
      cluster, split by grouping
      condition, downsampled.
      PNG format

  cmp_gr_clst_spl_cnd_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_clst_spl_cnd_res_*.pdf"
    doc: |
      Composition plot, colored by
      cluster, split by grouping
      condition, downsampled.
      PDF format

  cmp_gr_cnd_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_cnd_spl_clst_res_*.png"
    doc: |
      Composition plot, colored by
      grouping condition, split by
      cluster, downsampled.
      PNG format

  cmp_gr_cnd_spl_clst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_cnd_spl_clst_res_*.pdf"
    doc: |
      Composition plot, colored by
      grouping condition, split by
      cluster, downsampled.
      PDF format

  umap_spl_ph_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_ph_res_*.png"
    doc: |
      UMAP, colored by cluster,
      split by cell cycle phase.
      PNG format

  umap_spl_ph_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_ph_res_*.pdf"
    doc: |
      UMAP, colored by cluster,
      split by cell cycle phase.
      PDF format

  cmp_gr_ph_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_idnt.png"
    doc: |
      Composition plot, colored by
      cell cycle phase, split by
      dataset, downsampled.
      PNG format

  cmp_gr_ph_spl_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_idnt.pdf"
    doc: |
      Composition plot, colored by
      cell cycle phase, split by
      dataset, downsampled.
      PDF format

  cmp_gr_ph_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_ph_spl_clst_res_*.png"
    doc: |
      Composition plot, colored by
      cell cycle phase, split by
      cluster, downsampled.
      PNG format

  cmp_gr_ph_spl_clst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_ph_spl_clst_res_*.pdf"
    doc: |
      Composition plot, colored by
      cell cycle phase, split by
      cluster, downsampled.
      PDF format

  xpr_avg_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_avg_res_*.png"
    doc: |
      Gene expression dot plot.
      PNG format

  xpr_avg_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_avg_res_*.pdf"
    doc: |
      Gene expression dot plot.
      PDF format

  xpr_per_cell_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_[!sgnl_]*.png"
    doc: |
      UMAP, gene expression.
      PNG format

  xpr_per_cell_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_[!sgnl_]*.pdf"
    doc: |
      UMAP, gene expression.
      PDF format

  xpr_per_cell_sgnl_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_sgnl_*.png"
    doc: |
      UMAP, gene expression density.
      PNG format

  xpr_per_cell_sgnl_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_sgnl_*.pdf"
    doc: |
      UMAP, gene expression density.
      PDF format

  xpr_dnst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_dnst_res_*.png"
    doc: |
      Gene expression violin plot.
      PNG format

  xpr_dnst_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_dnst_res_*.pdf"
    doc: |
      Gene expression violin plot.
      PDF format

  cvrg_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cvrg_res_*.png"
    doc: |
      Fragments coverage.
      PNG format

  cvrg_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cvrg_res_*.pdf"
    doc: |
      Fragments coverage.
      PDF format

  xpr_htmp_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_htmp_res_*.png"
    doc: |
      Gene expression heatmap.
      PNG format

  xpr_htmp_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_htmp_res_*.pdf"
    doc: |
      Gene expression heatmap.
      PDF format

  xpr_htmp_res_tsv:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_htmp_res_*.tsv"
    doc: |
      Gene markers used for gene
      expression heatmap.
      TSV format

  gene_markers_tsv:
    type: File?
    outputBinding:
      glob: "*_gene_markers.tsv"
    doc: |
      Gene markers per cluster for
      all resolutions.
      TSV format

  peak_markers_tsv:
    type: File?
    outputBinding:
      glob: "*_peak_markers.tsv"
    doc: |
      Peak markers per cluster for
      all resolutions.
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


baseCommand: ["sc_wnn_cluster.R"]

stdout: sc_wnn_cluster_stdout.log
stderr: sc_wnn_cluster_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell WNN Cluster Analysis"
s:name: "Single-cell WNN Cluster Analysis"
s:alternateName: "Clusters multiome ATAC and RNA-Seq datasets, identifies gene markers and differentially accessible peaks"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-wnn-cluster.cwl
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
  Single-cell WNN Cluster Analysis

  Clusters multiome ATAC and RNA-Seq datasets, identifies gene markers
  and differentially accessible peaks.


s:about: |
  usage: sc_wnn_cluster.R
        [-h] --query QUERY
        [--rnadimensions [RNADIMENSIONS [RNADIMENSIONS ...]]]
        [--atacdimensions [ATACDIMENSIONS [ATACDIMENSIONS ...]]]
        [--algorithm {louvain,mult-louvain,slm,leiden}] [--uspread USPREAD]
        [--umindist UMINDIST] [--uneighbors UNEIGHBORS]
        [--umetric {euclidean,manhattan,chebyshev,minkowski,canberra,braycurtis,mahalanobis,wminkowski,seuclidean,cosine,correlation,haversine,hamming,jaccard,dice,russelrao,kulsinski,ll_dirichlet,hellinger,rogerstanimoto,sokalmichener,sokalsneath,yule}]
        [--umethod {uwot,uwot-learn,umap-learn}]
        [--resolution [RESOLUTION [RESOLUTION ...]]] [--fragments FRAGMENTS]
        [--genes [GENES [GENES ...]]] [--diffgenes] [--diffpeaks]
        [--rnalogfc RNALOGFC] [--rnaminpct RNAMINPCT] [--rnaonlypos]
        [--rnatestuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}]
        [--ataclogfc ATACLOGFC] [--atacminpct ATACMINPCT]
        [--atactestuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}]
        [--pdf] [--verbose] [--h5seurat] [--h5ad] [--cbbuild] [--output OUTPUT]
        [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
        [--cpus CPUS] [--memory MEMORY]

  Single-cell WNN Cluster Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include genes expression and chromatin
                          accessibility information stored in the RNA and ATAC
                          assays correspondingly. Additionally, 'pca',
                          'rnaumap', 'atac_lsi' and 'atacumap' dimensionality
                          reductions should be present.
    --rnadimensions [RNADIMENSIONS [RNADIMENSIONS ...]]
                          Dimensionality from the 'pca' reduction to use when
                          constructing weighted nearest-neighbor graph before
                          clustering (from 1 to 50). If single value N is
                          provided, use from 1 to N dimensions. If multiple
                          values are provided, subset to only selected
                          dimensions. Default: from 1 to 10
    --atacdimensions [ATACDIMENSIONS [ATACDIMENSIONS ...]]
                          Dimensionality from the 'atac_lsi' reduction to use
                          when constructing weighted nearest-neighbor graph
                          before clustering (from 1 to 50). If single value N is
                          provided, use from 2 to N dimensions. If multiple
                          values are provided, subset to only selected
                          dimensions. Default: from 2 to 10
    --algorithm {louvain,mult-louvain,slm,leiden}
                          Algorithm for modularity optimization when running
                          clustering. Default: louvain
    --uspread USPREAD     The effective scale of embedded points on UMAP. In
                          combination with '--mindist' it determines how
                          clustered/clumped the embedded points are. Default: 1
    --umindist UMINDIST   Controls how tightly the embedding is allowed compress
                          points together on UMAP. Larger values ensure embedded
                          points are moreevenly distributed, while smaller
                          values allow the algorithm to optimise more accurately
                          with regard to local structure. Sensible values are in
                          the range 0.001 to 0.5. Default: 0.3
    --uneighbors UNEIGHBORS
                          Determines the number of neighboring points used in
                          UMAP. Larger values will result in more global
                          structure being preserved at the loss of detailed
                          local structure. In general this parameter should
                          often be in the range 5 to 50. Default: 30
    --umetric {euclidean,manhattan,chebyshev,minkowski,canberra,braycurtis,mahalanobis,wminkowski,seuclidean,cosine,correlation,haversine,hamming,jaccard,dice,russelrao,kulsinski,ll_dirichlet,hellinger,rogerstanimoto,sokalmichener,sokalsneath,yule}
                          The metric to use to compute distances in high
                          dimensional space for UMAP. Default: cosine
    --umethod {uwot,uwot-learn,umap-learn}
                          UMAP implementation to run. If set to 'umap-learn' use
                          --umetric 'correlation' Default: uwot
    --resolution [RESOLUTION [RESOLUTION ...]]
                          Clustering resolution applied to the constructed
                          weighted nearest-neighbor graph. Can be set as an
                          array but only the first item from the list will be
                          used for cluster labels and gene/peak markers in the
                          UCSC Cell Browser when running with --cbbuild and
                          --diffgenes/--diffpeaks parameters. Default: 0.3, 0.5,
                          1.0
    --fragments FRAGMENTS
                          Count and barcode information for every ATAC fragment
                          used in the loaded Seurat object. File should be saved
                          in TSV format with tbi-index file.
    --genes [GENES [GENES ...]]
                          Genes of interest to build gene expression and Tn5
                          insertion frequency plots for the nearest peaks. If '
                          --fragments' is not provided only gene expression
                          plots will be built. Default: None
    --diffgenes           Identify differentially expressed genes (putative gene
                          markers) between each pair of clusters for all
                          resolutions. Default: false
    --diffpeaks           Identify differentially accessible peaks between each
                          pair of clusters for all resolutions. Default: false
    --rnalogfc RNALOGFC   For putative gene markers identification include only
                          those genes that on average have log fold change
                          difference in expression between every tested pair of
                          clusters not lower than this value. Ignored if '--
                          diffgenes' is not set. Default: 0.25
    --rnaminpct RNAMINPCT
                          For putative gene markers identification include only
                          those genes that are detected in not lower than this
                          fraction of cells in either of the two tested
                          clusters. Ignored if '--diffgenes' is not set.
                          Default: 0.1
    --rnaonlypos          For putative gene markers identification return only
                          positive markers. Ignored if '--diffgenes' is not set.
                          Default: false
    --rnatestuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}
                          Statistical test to use for putative gene markers
                          identification. Ignored if '--diffgenes' is not set.
                          Default: wilcox
    --ataclogfc ATACLOGFC
                          For differentially accessible peaks identification
                          include only those peaks that on average have log fold
                          change difference in the chromatin accessibility
                          between every tested pair of clusters not lower than
                          this value. Ignored if '--diffpeaks' is not set.
                          Default: 0.25
    --atacminpct ATACMINPCT
                          For differentially accessible peaks identification
                          include only those peaks that are detected in not
                          lower than this fraction of cells in either of the two
                          tested clusters. Ignored if '--diffpeaks' is not set.
                          Default: 0.05
    --atactestuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}
                          Statistical test to use for differentially accessible
                          peaks identification. Ignored if '--diffpeaks' is not
                          set. Default: LR
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