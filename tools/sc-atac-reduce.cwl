cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.42


inputs:

  query_data_rds:
    type: File
    inputBinding:
      prefix: "--query"
    doc: |
      Path to the RDS file to load Seurat object from. This file should include
      chromatin accessibility information stored in the ATAC assay.

  datasets_metadata:
    type: File?
    inputBinding:
      prefix: "--metadata"
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat
      object metadata with categorical values using samples
      identities. First column - 'library_id' should
      correspond to all unique values from the 'new.ident'
      column of the loaded Seurat object. If any of the
      provided in this file columns are already present in
      the Seurat object metadata, they will be overwritten.
      When combined with --barcodes parameter, first the
      metadata will be extended, then barcode filtering will
      be applied.
      Default: no extra metadata is added

  barcodes_data:
    type: File?
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the TSV/CSV file to optionally prefilter and
      extend Seurat object metadata be selected barcodes.
      First column should be named as 'barcode'. If file
      includes any other columns they will be added to the
      Seurat object metadata ovewriting the existing ones if
      those are present.
      Default: all cells used, no extra metadata is added

  normalization_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "log-tfidf"
      - "tf-logidf"
      - "logtf-logidf"
      - "idf"
    inputBinding:
      prefix: "--norm"
    doc: |
      TF-IDF normalization method applied to chromatin
      accessibility counts. log-tfidf - Stuart & Butler et
      al. 2019, tf-logidf - Cusanovich & Hill et al. 2018,
      logtf-logidf - Andrew Hill, idf - 10x Genomics,
      Default: log-tfidf

  integration_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "signac"
      - "harmony"
      - "none"
    inputBinding:
      prefix: "--ntgr"
    doc: |
      Integration method used for joint analysis of multiple
      datasets.
      Default: signac

  integrate_by:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--ntgrby"
    doc: |
      Column(s) from the Seurat object metadata to define
      the variable(s) that should be integrated out when
      running multiple datasets integration with harmony.
      May include columns from the extra metadata added with
      --metadata parameter. Ignored if --ntgr is not set to
      harmony.
      Default: new.ident

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
    type: int?
    inputBinding:
      prefix: "--dimensions"
    doc: |
      Dimensionality to use for datasets integration (if provided RDS file includes
      multiple datasets and --ntgr is not set to 'none') and UMAP projection.
      (from 2 to 50). First LSI component is always excluded.
      Default: 10

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

  qc_dim_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_qc_dim_corr.png"
    doc: |
      Correlation between QC metrics and LSI components.
      PNG format.

  qc_dim_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_qc_dim_corr.pdf"
    doc: |
      Correlation between QC metrics and LSI components.
      PDF format.

  umap_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_qc_mtrcs.png"
    doc: |
      UMAP, QC metrics.
      PNG format.

  umap_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_qc_mtrcs.pdf"
    doc: |
      UMAP, QC metrics.
      PDF format.

  umap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap.png"
    doc: |
      UMAP, colored by dataset.
      PNG format.

  umap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap.pdf"
    doc: |
      UMAP, colored by dataset.
      PDF format.

  umap_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt.png"
    doc: |
      UMAP, split by dataset.
      PNG format.

  umap_spl_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt.pdf"
    doc: |
      UMAP, split by dataset.
      PDF format.

  umap_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_cnd.png"
    doc: |
      UMAP, colored by dataset, split
      by grouping condition.
      PNG format.

  umap_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_cnd.pdf"
    doc: |
      UMAP, colored by dataset, split
      by grouping condition.
      PDF format.

  umap_spl_frgm_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_frgm.png"
    doc: |
      UMAP, colored by dataset, split
      by ATAC fragments in peaks per cell.
      PNG format.

  umap_spl_frgm_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_frgm.pdf"
    doc: |
      UMAP, colored by dataset, split
      by ATAC fragments in peaks per cell.
      PDF format.

  umap_spl_peak_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_peak.png"
    doc: |
      UMAP, colored by dataset, split
      by peaks per cell.
      PNG format.

  umap_spl_peak_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_peak.pdf"
    doc: |
      UMAP, colored by dataset, split
      by peaks per cell.
      PDF format.

  umap_spl_tss_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_tss.png"
    doc: |
      UMAP, colored by dataset, split
      by TSS enrichment score.
      PNG format.

  umap_spl_tss_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_tss.pdf"
    doc: |
      UMAP, colored by dataset, split
      by TSS enrichment score.
      PDF format.

  umap_spl_ncls_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_ncls.png"
    doc: |
      UMAP, colored by dataset, split
      by nucleosome signal.
      PNG format.

  umap_spl_ncls_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_ncls.pdf"
    doc: |
      UMAP, colored by dataset, split
      by nucleosome signal.
      PDF format.

  umap_spl_frip_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_frip.png"
    doc: |
      UMAP, colored by dataset,
      split by FRiP.
      PNG format.

  umap_spl_frip_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_frip.pdf"
    doc: |
      UMAP, colored by dataset,
      split by FRiP.
      PDF format.

  umap_spl_blck_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_blck.png"
    doc: |
      UMAP, colored by dataset, split
      by blacklist fraction.
      PNG format.

  umap_spl_blck_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_blck.pdf"
    doc: |
      UMAP, colored by dataset, split
      by blacklist fraction.
      PDF format.

  umap_gr_cnd_spl_frgm_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_cnd_spl_frgm.png"
    doc: |
      UMAP, colored by grouping condition,
      split by ATAC fragments in peaks per cell.
      PNG format.

  umap_gr_cnd_spl_frgm_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_gr_cnd_spl_frgm.pdf"
    doc: |
      UMAP, colored by grouping condition,
      split by ATAC fragments in peaks per cell.
      PDF format.

  umap_gr_cnd_spl_peak_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_cnd_spl_peak.png"
    doc: |
      UMAP, colored by grouping condition,
      split by peaks per cell.
      PNG format.

  umap_gr_cnd_spl_peak_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_gr_cnd_spl_peak.pdf"
    doc: |
      UMAP, colored by grouping condition,
      split by peaks per cell.
      PDF format.

  umap_gr_cnd_spl_tss_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_cnd_spl_tss.png"
    doc: |
      UMAP, colored by grouping condition,
      split by TSS enrichment score.
      PNG format.

  umap_gr_cnd_spl_tss_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_gr_cnd_spl_tss.pdf"
    doc: |
      UMAP, colored by grouping condition,
      split by TSS enrichment score.
      PDF format.

  umap_gr_cnd_spl_ncls_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_cnd_spl_ncls.png"
    doc: |
      UMAP, colored by grouping condition,
      split by nucleosome signal.
      PNG format.

  umap_gr_cnd_spl_ncls_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_gr_cnd_spl_ncls.pdf"
    doc: |
      UMAP, colored by grouping condition,
      split by nucleosome signal.
      PDF format.

  umap_gr_cnd_spl_frip_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_cnd_spl_frip.png"
    doc: |
      UMAP, colored by grouping condition,
      split by FRiP.
      PNG format.

  umap_gr_cnd_spl_frip_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_gr_cnd_spl_frip.pdf"
    doc: |
      UMAP, colored by grouping condition,
      split by FRiP.
      PDF format.

  umap_gr_cnd_spl_blck_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_cnd_spl_blck.png"
    doc: |
      UMAP, colored by grouping condition,
      split by blacklist fraction.
      PNG format.

  umap_gr_cnd_spl_blck_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_gr_cnd_spl_blck.pdf"
    doc: |
      UMAP, colored by grouping condition,
      split by blacklist fraction.
      PDF format.

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
- valueFrom: $(inputs.export_html_report?["/usr/local/bin/sc_report_wrapper.R", "/usr/local/bin/sc_atac_reduce.R"]:"/usr/local/bin/sc_atac_reduce.R")


stdout: sc_atac_reduce_stdout.log
stderr: sc_atac_reduce_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-Cell ATAC-Seq Dimensionality Reduction Analysis"
s:name: "Single-Cell ATAC-Seq Dimensionality Reduction Analysis"
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
  Single-Cell ATAC-Seq Dimensionality Reduction Analysis

  Integrates multiple single-cell ATAC-Seq datasets, reduces dimensionality using LSI.


s:about: |
  usage: /usr/local/bin/sc_atac_reduce.R [-h] --query QUERY
                                        [--metadata METADATA]
                                        [--barcodes BARCODES]
                                        [--norm {log-tfidf,tf-logidf,logtf-logidf,idf}]
                                        [--ntgr {signac,harmony,none}]
                                        [--ntgrby [NTGRBY [NTGRBY ...]]]
                                        [--minvarpeaks MINVARPEAKS]
                                        [--dimensions DIMENSIONS]
                                        [--uspread USPREAD]
                                        [--umindist UMINDIST]
                                        [--uneighbors UNEIGHBORS]
                                        [--umetric {euclidean,manhattan,chebyshev,minkowski,canberra,braycurtis,mahalanobis,wminkowski,seuclidean,cosine,correlation,haversine,hamming,jaccard,dice,russelrao,kulsinski,ll_dirichlet,hellinger,rogerstanimoto,sokalmichener,sokalsneath,yule}]
                                        [--umethod {uwot,uwot-learn,umap-learn}]
                                        [--pdf] [--verbose] [--h5seurat]
                                        [--h5ad] [--cbbuild] [--output OUTPUT]
                                        [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
                                        [--cpus CPUS] [--memory MEMORY]
                                        [--seed SEED]

  Single-Cell ATAC-Seq Dimensionality Reduction Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include chromatin accessibility
                          information stored in the ATAC assay.
    --metadata METADATA   Path to the TSV/CSV file to optionally extend Seurat
                          object metadata with categorical values using samples
                          identities. First column - 'library_id' should
                          correspond to all unique values from the 'new.ident'
                          column of the loaded Seurat object. If any of the
                          provided in this file columns are already present in
                          the Seurat object metadata, they will be overwritten.
                          When combined with --barcodes parameter, first the
                          metadata will be extended, then barcode filtering will
                          be applied. Default: no extra metadata is added
    --barcodes BARCODES   Path to the TSV/CSV file to optionally prefilter and
                          extend Seurat object metadata be selected barcodes.
                          First column should be named as 'barcode'. If file
                          includes any other columns they will be added to the
                          Seurat object metadata ovewriting the existing ones if
                          those are present. Default: all cells used, no extra
                          metadata is added
    --norm {log-tfidf,tf-logidf,logtf-logidf,idf}
                          TF-IDF normalization method applied to chromatin
                          accessibility counts. log-tfidf - Stuart & Butler et
                          al. 2019, tf-logidf - Cusanovich & Hill et al. 2018,
                          logtf-logidf - Andrew Hill, idf - 10x Genomics,
                          Default: log-tfidf
    --ntgr {signac,harmony,none}
                          Integration method used for joint analysis of multiple
                          datasets. Default: signac
    --ntgrby [NTGRBY [NTGRBY ...]]
                          Column(s) from the Seurat object metadata to define
                          the variable(s) that should be integrated out when
                          running multiple datasets integration with harmony.
                          May include columns from the extra metadata added with
                          --metadata parameter. Ignored if --ntgr is not set to
                          harmony. Default: new.ident
    --minvarpeaks MINVARPEAKS
                          Minimum percentile for identifying the top most common
                          peaks as highly variable. For example, setting to 5
                          will use the the top 95 percent most common among all
                          cells peaks as highly variable. These peaks are used
                          for datasets integration, scaling and dimensionality
                          reduction. Default: 0 (use all available peaks)
    --dimensions DIMENSIONS
                          Dimensionality to use for datasets integration (if
                          provided RDS file includes multiple datasets and
                          --ntgr is not set to 'none') and UMAP projection.
                          (from 2 to 50). First LSI component is always
                          excluded. Default: 10
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