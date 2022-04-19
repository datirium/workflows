cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $(inputs.vector_memory_limit * 1000000000)


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.4


inputs:

  query_data_rds:
    type: File
    inputBinding:
      prefix: "--query"
    doc: |
      Path to the RDS file to load Seurat object from. This file should include
      chromatin accessibility information stored in the ATAC assay.

  barcodes_data:
    type: File?
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the headerless TSV/CSV file with the list of barcodes to select
      cells of interest (one barcode per line). Prefilters loaded Seurat object
      to include only specific set of cells.
      Default: use all cells.

  integration_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "signac"
      - "none"
    inputBinding:
      prefix: "--ntgr"
    doc: |
      Integration method used for joint analysis of multiple datasets. Automatically
      set to 'none' if loaded Suerat object includes only one dataset.
      Default: signac

  minimum_var_peaks_perc:
    type: int?
    inputBinding:
      prefix: "--minvarpeaks"
    doc: |
      Minimum percentile for identifying the top most common peaks as highly variable.
      For example, setting to 5 will use the the top 95 percent most common among all cells
      peaks as highly variable. These peaks are used for datasets integration, scaling
      and dimensionality reduction.
      Default: 0 (use all available peaks)

  dimensions:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--dimensions"
    doc: |
      Dimensionality to use for datasets integration and UMAP projection (from 2 to 50).
      If single value N is provided, use from 2 to N LSI components. If multiple values are
      provided, subset to only selected LSI components.
      Default: from 2 to 10

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

  export_pdf_plots:
    type: boolean?
    inputBinding:
      prefix: "--pdf"
    doc: |
      Export plots in PDF.
      Default: false

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

  qc_dim_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_qc_dim_corr.png"
    doc: |
      Correlation plots between QC metrics and cells LSI dimensions.
      PNG format

  qc_dim_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_qc_dim_corr.pdf"
    doc: |
      Correlation plots between QC metrics and cells LSI dimensions.
      PDF format

  umap_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_qc_mtrcs.png"
    doc: |
      QC metrics on cells UMAP.
      PNG format

  umap_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_qc_mtrcs.pdf"
    doc: |
      QC metrics on cells UMAP.
      PDF format

  umap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap.png"
    doc: |
      Cells UMAP.
      PNG format

  umap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap.pdf"
    doc: |
      Cells UMAP.
      PDF format

  umap_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt.png"
    doc: |
      Split by dataset cells UMAP.
      PNG format

  umap_spl_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt.pdf"
    doc: |
      Split by dataset cells UMAP.
      PDF format

  umap_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_cnd.png"
    doc: |
      Split by grouping condition cells UMAP.
      PNG format

  umap_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_cnd.pdf"
    doc: |
      Split by grouping condition cells UMAP.
      PDF format

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

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sc_atac_reduce.R"]

stdout: sc_atac_reduce_stdout.log
stderr: sc_atac_reduce_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell ATAC-Seq Dimensionality Reduction Analysis"
s:name: "Single-cell ATAC-Seq Dimensionality Reduction Analysis"
s:alternateName: "Integrates multiple single-cell ATAC-Seq datasets, reduces dimensionality using LSI"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-atac-reduce.cwl
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
  Single-cell ATAC-Seq Dimensionality Reduction Analysis
  ====================================================================================
  Integrates multiple single-cell ATAC-Seq datasets, reduces dimensionality using LSI.


s:about: |
  usage: sc_atac_reduce.R
        [-h] --query QUERY [--barcodes BARCODES] [--ntgr {signac,none}]
        [--minvarpeaks MINVARPEAKS] [--dimensions [DIMENSIONS ...]]
        [--uspread USPREAD] [--umindist UMINDIST] [--uneighbors UNEIGHBORS]
        [--umetric {euclidean,manhattan,chebyshev,minkowski,canberra,braycurtis,mahalanobis,wminkowski,seuclidean,cosine,correlation,haversine,hamming,jaccard,dice,russelrao,kulsinski,ll_dirichlet,hellinger,rogerstanimoto,sokalmichener,sokalsneath,yule}]
        [--umethod {uwot,uwot-learn,umap-learn}] [--pdf] [--verbose]
        [--h5seurat] [--cbbuild] [--output OUTPUT] [--cpus CPUS]
        [--memory MEMORY]

  Single-cell ATAC-Seq Dimensionality Reduction Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include chromatin accessibility
                          information stored in the ATAC assay.
    --barcodes BARCODES   Path to the headerless TSV/CSV file with the list of
                          barcodes to select cells of interest (one barcode per
                          line). Prefilters loaded Seurat object to include only
                          specific set of cells. Default: use all cells.
    --ntgr {signac,none}  Integration method used for joint analysis of multiple
                          datasets. Automatically set to 'none' if loaded Suerat
                          object includes only one dataset. Default: signac
    --minvarpeaks MINVARPEAKS
                          Minimum percentile for identifying the top most common
                          peaks as highly variable. For example, setting to 5
                          will use the the top 95 percent most common among all
                          cells peaks as highly variable. These peaks are used
                          for datasets integration, scaling and dimensionality
                          reduction. Default: 0 (use all available peaks)
    --dimensions [DIMENSIONS ...]
                          Dimensionality to use for datasets integration and
                          UMAP projection (from 2 to 50). If single value N is
                          provided, use from 2 to N LSI components. If multiple
                          values are provided, subset to only selected LSI
                          components. Default: from 2 to 10
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
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --h5seurat            Save Seurat data to h5seurat file. Default: false
    --cbbuild             Export results to UCSC Cell Browser. Default: false
    --output OUTPUT       Output prefix. Default: ./sc
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32