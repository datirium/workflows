cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.41


inputs:

  query_data_rds:
    type: File
    inputBinding:
      prefix: "--query"
    doc: |
      Path to the RDS file to load Seurat object from. This file should include
      chromatin accessibility information stored in the ATAC assay, as well as
      'atac_lsi' and 'atacumap' dimensionality reductions applied to that assay.

  dimensions:
    type: int?
    inputBinding:
      prefix: "--dimensions"
    doc: |
      Dimensionality to use when constructing nearest-neighbor graph before clustering
      (from 2 to 50). First LSI component is always excluded unless the provided RDS
      file consists of multiple datasets integrated with Harmony.
      Default: 10

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
      Default: slm

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
      for cluster labels and peak markers in the UCSC Cell Browser when running
      with --cbbuild and --diffpeaks parameters.
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
      Genes of interest to build Tn5 insertion frequency plots for the nearest peaks.
      If loaded Seurat object includes genes expression information in the RNA assay
      it will be additionally shown on the right side of the plots.
      Ignored if '--fragments' is not provided.
      Default: None

  cvrg_upstream_bp:
    type: int?
    inputBinding:
      prefix: "--upstream"
    doc: |
      Number of bases to extend the genome coverage region for
      a specific gene upstream. Ignored if --genes or --fragments
      parameters are not provided. Default: 2500

  cvrg_downstream_bp:
    type: int?
    inputBinding:
      prefix: "--downstream"
    doc: |
      Number of bases to extend the genome coverage region for
      a specific gene downstream. Ignored if --genes or --fragments
      parameters are not provided. Default: 2500

  identify_diff_peaks:
    type: boolean?
    inputBinding:
      prefix: "--diffpeaks"
    doc: |
      Identify differentially accessible peaks between each pair of clusters for all resolutions.
      Default: false

  minimum_logfc:
    type: float?
    inputBinding:
      prefix: "--logfc"
    doc: |
      For differentially accessible peaks identification include only those peaks that
      on average have log fold change difference in the chromatin accessibility between
      every tested pair of clusters not lower than this value. Ignored if '--diffpeaks'
      is not set.
      Default: 0.25

  minimum_pct:
    type: float?
    inputBinding:
      prefix: "--minpct"
    doc: |
      For differentially accessible peaks identification include only those peaks that
      are detected in not lower than this fraction of cells in either of the two tested
      clusters. Ignored if '--diffpeaks' is not set.
      Default: 0.05

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
      Save raw counts from the ATAC assay to h5ad file.
      Default: false

  export_ucsc_cb:
    type: boolean?
    inputBinding:
      prefix: "--cbbuild"
    doc: |
      Export results to UCSC Cell Browser. Default: false

  export_html_report:
    type: boolean?
    default: false
    doc: |
      Export tehcnical report. HTML format.
      Note, stdout will be less informative.
      Default: false

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

  seed:
    type: int?
    inputBinding:
      prefix: "--seed"
    doc: |
      Seed number for random values.
      Default: 42


outputs:

  cell_cnts_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cell_cnts_gr_clst_res_*.png"
    doc: |
      Number of cells per cluster.
      All cells; all resolutions.
      PNG format.

  tss_frgm_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_tss_frgm_spl_clst_res_*.png"
    doc: |
      TSS enrichment score vs ATAC
      fragments in peaks per cell.
      Split by cluster; all cells;
      all resolutions.
      PNG format.

  atacdbl_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_atacdbl_gr_clst_res_*.png"
    doc: |
      Percentage of ATAC doublets per cluster.
      All cells; all resolutions.
      PNG format.

  qc_mtrcs_dnst_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_qc_mtrcs_dnst_gr_clst_res_*.png"
    doc: |
      Distribution of QC metrics per cell
      colored by cluster.
      All cells; all resolutions.
      PNG format.

  umap_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_gr_clst_res_*.png"
    doc: |
      UMAP colored by cluster.
      All cells; all resolutions.
      PNG format.

  slh_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_slh_gr_clst_res_*.png"
    doc: |
      Silhouette scores.
      All cells; all resolutions.
      PNG format.

  umap_gr_clst_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_gr_clst_spl_idnt_res_*.png"
    doc: |
      UMAP colored by cluster.
      Split by dataset; downsampled to the
      smallest dataset; all resolutions.
      PNG format.

  cmp_gr_clst_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_clst_spl_idnt_res_*.png"
    doc: |
      Composition plot colored by cluster.
      Split by dataset; downsampled to the
      smallest dataset; all resolutions.
      PNG format.

  cmp_gr_idnt_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_idnt_spl_clst_res_*.png"
    doc: |
      Composition plot colored by dataset.
      Split by cluster; downsampled to the
      smallest dataset; all resolutions.
      PNG format.

  umap_gr_clst_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_gr_clst_spl_cnd_res_*.png"
    doc: |
      UMAP colored by cluster.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group; all resolutions.
      PNG format.

  cmp_gr_clst_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_clst_spl_cnd_res_*.png"
    doc: |
      Composition plot colored by cluster.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group; all resolutions.
      PNG format.

  cmp_gr_cnd_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cmp_gr_cnd_spl_clst_res_*.png"
    doc: |
      Composition plot colored by grouping condition.
      Split by cluster; first downsampled to the
      smallest dataset, then downsampled to the
      smallest group; all resolutions.
      PNG format.

  cvrg_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cvrg_res_*.png"
    doc: |
      ATAC fragment coverage.
      All genes of interest; all resolutions.
      PNG format.

  all_plots_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*.pdf"
    doc: |
      All generated plots.
      PDF format.

  peak_markers_tsv:
    type: File?
    outputBinding:
      glob: "*_peak_markers.tsv"
    doc: |
      Peak markers.
      All resolutions.
      TSV format.

  ucsc_cb_config_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser"
    doc: |
      UCSC Cell Browser configuration data.

  ucsc_cb_html_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser/html_data"
    doc: |
      UCSC Cell Browser html data.

  ucsc_cb_html_file:
    type: File?
    outputBinding:
      glob: "*_cellbrowser/html_data/index.html"
    doc: |
      UCSC Cell Browser html index.

  seurat_data_rds:
    type: File
    outputBinding:
      glob: "*_data.rds"
    doc: |
      Seurat object.
      RDS format.

  seurat_data_h5seurat:
    type: File?
    outputBinding:
      glob: "*_data.h5seurat"
    doc: |
      Seurat object.
      h5Seurat format.

  seurat_data_h5ad:
    type: File?
    outputBinding:
      glob: "*_counts.h5ad"
    doc: |
      Seurat object.
      H5AD format.

  sc_report_html_file:
    type: File?
    outputBinding:
      glob: "sc_report.html"
    doc: |
      Tehcnical report.
      HTML format.

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["Rscript"]
arguments:
- valueFrom: $(inputs.export_html_report?["/usr/local/bin/sc_report_wrapper.R", "/usr/local/bin/sc_atac_cluster.R"]:"/usr/local/bin/sc_atac_cluster.R")


stdout: sc_atac_cluster_stdout.log
stderr: sc_atac_cluster_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-Cell ATAC-Seq Cluster Analysis"
s:name: "Single-Cell ATAC-Seq Cluster Analysis"
s:alternateName: "Clusters single-cell ATAC-Seq datasets, identifies differentially accessible peaks"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-atac-cluster.cwl
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
  Single-Cell ATAC-Seq Cluster Analysis

  Clusters single-cell ATAC-Seq datasets, identifies differentially
  accessible peaks.


s:about: |
  usage: /usr/local/bin/sc_atac_cluster.R [-h] --query QUERY
                                          [--dimensions DIMENSIONS]
                                          [--ametric {euclidean,cosine,manhattan,hamming}]
                                          [--algorithm {louvain,mult-louvain,slm,leiden}]
                                          [--resolution [RESOLUTION [RESOLUTION ...]]]
                                          [--fragments FRAGMENTS]
                                          [--genes [GENES [GENES ...]]]
                                          [--upstream UPSTREAM]
                                          [--downstream DOWNSTREAM]
                                          [--diffpeaks] [--logfc LOGFC]
                                          [--minpct MINPCT]
                                          [--testuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}]
                                          [--pdf] [--verbose] [--h5seurat]
                                          [--h5ad] [--cbbuild] [--output OUTPUT]
                                          [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
                                          [--cpus CPUS] [--memory MEMORY]
                                          [--seed SEED]

  Single-Cell ATAC-Seq Cluster Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include chromatin accessibility
                          information stored in the ATAC assay, as well as
                          'atac_lsi' and 'atacumap' dimensionality reductions
                          applied to that assay.
    --dimensions DIMENSIONS
                          Dimensionality to use when constructing nearest-
                          neighbor graph before clustering (from 2 to 50). First
                          LSI component is always excluded unless the provided
                          RDS file consists of multiple datasets integrated with
                          Harmony. Default: 10
    --ametric {euclidean,cosine,manhattan,hamming}
                          Distance metric used when constructing nearest-
                          neighbor graph before clustering. Default: euclidean
    --algorithm {louvain,mult-louvain,slm,leiden}
                          Algorithm for modularity optimization when running
                          clustering. Default: slm
    --resolution [RESOLUTION [RESOLUTION ...]]
                          Clustering resolution applied to the constructed
                          nearest-neighbor graph. Can be set as an array but
                          only the first item from the list will be used for
                          cluster labels and peak markers in the UCSC Cell
                          Browser when running with --cbbuild and --diffpeaks
                          parameters. Default: 0.3, 0.5, 1.0
    --fragments FRAGMENTS
                          Count and barcode information for every ATAC fragment
                          used in the loaded Seurat object. File should be saved
                          in TSV format with tbi-index file.
    --genes [GENES [GENES ...]]
                          Genes of interest to build Tn5 insertion frequency
                          plots for the nearest peaks. If loaded Seurat object
                          includes genes expression information in the RNA assay
                          it will be additionally shown on the right side of the
                          plots. Ignored if '--fragments' is not provided.
                          Default: None
    --upstream UPSTREAM   Number of bases to extend the genome coverage region
                          for a specific gene upstream. Ignored if --genes or
                          --fragments parameters are not provided. Default: 2500
    --downstream DOWNSTREAM
                          Number of bases to extend the genome coverage region
                          for a specific gene downstream. Ignored if --genes or
                          --fragments parameters are not provided. Default: 2500
    --diffpeaks           Identify differentially accessible peaks between each
                          pair of clusters for all resolutions. Default: false
    --logfc LOGFC         For differentially accessible peaks identification
                          include only those peaks that on average have log fold
                          change difference in the chromatin accessibility
                          between every tested pair of clusters not lower than
                          this value. Ignored if '--diffpeaks' is not set.
                          Default: 0.25
    --minpct MINPCT       For differentially accessible peaks identification
                          include only those peaks that are detected in not
                          lower than this fraction of cells in either of the two
                          tested clusters. Ignored if '--diffpeaks' is not set.
                          Default: 0.05
    --testuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}
                          Statistical test to use for differentially accessible
                          peaks identification. Ignored if '--diffpeaks' is not
                          set. Default: LR
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --h5seurat            Save Seurat data to h5seurat file. Default: false
    --h5ad                Save raw counts from the ATAC assay to h5ad file.
                          Default: false
    --cbbuild             Export results to UCSC Cell Browser. Default: false
    --output OUTPUT       Output prefix. Default: ./sc
    --theme {gray,bw,linedraw,light,dark,minimal,classic,void}
                          Color theme for all generated plots. Default: classic
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32
    --seed SEED           Seed number for random values. Default: 42