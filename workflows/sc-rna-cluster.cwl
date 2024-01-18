cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var split_features = function(line) {
          function get_unique(value, index, self) {
            return self.indexOf(value) === index && value != "";
          }
          let splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };
    - var split_numbers = function(line) {
          let splitted_line = line?line.split(/[\s,]+/).map(parseFloat):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };


'sd:upstream':
  sc_tools_sample:
  - "sc-rna-cluster.cwl"
  - "sc-atac-cluster.cwl"
  - "sc-rna-reduce.cwl"
  - "sc-atac-reduce.cwl"


inputs:

  alias:
    type: string
    label: "Analysis name"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Single-cell Analysis with PCA Transformed RNA-Seq Datasets"
    doc: |
      Analysis that includes single-cell
      multiome RNA and ATAC-Seq or just
      RNA-Seq datasets run through
      "Single-Cell RNA-Seq Dimensionality
      Reduction Analysis" at any of the
      processing stages.
    'sd:upstreamSource': "sc_tools_sample/seurat_data_rds"
    'sd:localLabel': true

  dimensions:
    type: int?
    default: 40
    label: "Target dimensionality"
    doc: |
      Number of principal components to be
      used in constructing nearest-neighbor
      graph as part of the clustering
      algorithm. Accepted values range from
      1 to 50.
      Default: 40

  resolution:
    type: float?
    default: 0.3
    label: "Clustering resolution"
    doc: |
      Resolution to define the "granularity"
      of the clustered data. Larger values
      lead to a bigger number of clusters.
      Optimal resolution often increases
      with the number of cells. For a dataset
      of 3K cells, the value within 0.4-1.2
      range usually returns good results.
      Default: 0.3

  identify_diff_genes:
    type: boolean?
    default: true
    label: "Find gene markers"
    doc: |
      Identify upregulated genes in each
      cluster compared to all other cells.
      Include only genes that are expressed
      in at least 10% of the cells coming
      from either current cluster or from
      all other clusters together.
      Exclude cells with log2FoldChange
      values less than 0.25. Use Wilcoxon
      Rank Sum test to calculate P-values.
      Keep only genes with P-values lower
      than 0.01. Adjust P-values for multiple
      comparisons using Bonferroni correction.
      Default: true

  genes_of_interest:
    type: string?
    default: null
    label: "Genes of interest"
    doc: |
      Comma or space separated list of
      genes of interest to visualize
      expression.
      Default: None

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
    default: "classic"
    label: "Plots color theme"
    doc: |
      Color theme for all plots saved
      as PNG files.
      Default: classic
    "sd:layout":
      advanced: true

  threads:
    type:
    - "null"
    - type: enum
      symbols:
      - "1"
      - "2"
      - "3"
      - "4"
      - "5"
      - "6"
    default: "1"
    label: "Cores/CPUs"
    doc: |
      Parallelization parameter to define the
      number of cores/CPUs that can be utilized
      simultaneously.
      Default: 1
    "sd:layout":
      advanced: true


outputs:

  umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/umap_res_plot_png
    label: "UMAP, colored by cluster"
    doc: |
      UMAP, colored by cluster
    'sd:visualPlugins':
    - image:
        tab: 'Per cluster'
        Caption: 'UMAP, colored by cluster'

  slh_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/slh_res_plot_png
    label: "Silhouette scores"
    doc: |
      Silhouette scores
    'sd:visualPlugins':
    - image:
        tab: 'Per cluster'
        Caption: 'Silhouette scores'

  umap_spl_ph_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/umap_spl_ph_res_plot_png
    label: "UMAP, colored by cluster, split by cell cycle phase"
    doc: |
      UMAP, colored by cluster,
      split by cell cycle phase
    'sd:visualPlugins':
    - image:
        tab: 'Per cluster'
        Caption: 'UMAP, colored by cluster, split by cell cycle phase'

  cmp_gr_ph_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/cmp_gr_ph_spl_clst_res_plot_png
    label: "Composition plot, colored by cell cycle phase, split by cluster, downsampled"
    doc: |
      Composition plot, colored by
      cell cycle phase, split by
      cluster, downsampled
    'sd:visualPlugins':
    - image:
        tab: 'Per cluster'
        Caption: 'Composition plot, colored by cell cycle phase, split by cluster, downsampled'

  umap_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/umap_spl_idnt_res_plot_png
    label: "UMAP, colored by cluster, split by dataset"
    doc: |
      UMAP, colored by cluster,
      split by dataset
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'UMAP, colored by cluster, split by dataset'

  cmp_gr_clst_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/cmp_gr_clst_spl_idnt_res_plot_png
    label: "Composition plot, colored by cluster, split by dataset, downsampled"
    doc: |
      Composition plot, colored by
      cluster, split by dataset,
      downsampled
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Composition plot, colored by cluster, split by dataset, downsampled'

  cmp_gr_idnt_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/cmp_gr_idnt_spl_clst_res_plot_png
    label: "Composition plot, colored by dataset, split by cluster, downsampled"
    doc: |
      Composition plot, colored by
      dataset, split by cluster,
      downsampled
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Composition plot, colored by dataset, split by cluster, downsampled'

  cmp_gr_ph_spl_idnt_plot_png:
    type: File?
    outputSource: sc_rna_cluster/cmp_gr_ph_spl_idnt_plot_png
    label: "Composition plot, colored by cell cycle phase, split by dataset, downsampled"
    doc: |
      Composition plot, colored by
      cell cycle phase, split by
      dataset, downsampled
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Composition plot, colored by cell cycle phase, split by dataset, downsampled'

  umap_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/umap_spl_cnd_res_plot_png
    label: "UMAP, colored by cluster, split by grouping condition"
    doc: |
      UMAP, colored by cluster,
      split by grouping condition
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'UMAP, colored by cluster, split by grouping condition'

  cmp_gr_clst_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/cmp_gr_clst_spl_cnd_res_plot_png
    label: "Composition plot, colored by cluster, split by grouping condition, downsampled"
    doc: |
      Composition plot, colored by
      cluster, split by grouping
      condition, downsampled
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Composition plot, colored by cluster, split by grouping condition, downsampled'

  cmp_gr_cnd_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/cmp_gr_cnd_spl_clst_res_plot_png
    label: "Composition plot, colored by grouping condition, split by cluster, downsampled"
    doc: |
      Composition plot, colored by
      grouping condition, split by
      cluster, downsampled
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Composition plot, colored by grouping condition, split by cluster, downsampled'

  xpr_avg_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/xpr_avg_res_plot_png
    label: "Gene expression dot plot"
    doc: |
      Gene expression dot plot
    'sd:visualPlugins':
    - image:
        tab: 'Genes of interest'
        Caption: 'Gene expression dot plot'

  xpr_dnst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/xpr_dnst_res_plot_png
    label: "Gene expression violin plot"
    doc: |
      Gene expression violin plot
    'sd:visualPlugins':
    - image:
        tab: 'Genes of interest'
        Caption: 'Gene expression violin plot'

  xpr_per_cell_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/xpr_per_cell_plot_png
    label: "UMAP, gene expression"
    doc: |
      UMAP, gene expression
    'sd:visualPlugins':
    - image:
        tab: 'Genes of interest'
        Caption: 'UMAP, gene expression'

  xpr_htmp_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/xpr_htmp_res_plot_png
    label: "Gene expression heatmap"
    doc: |
      Gene expression heatmap
    'sd:visualPlugins':
    - image:
        tab: 'Heatmap'
        Caption: 'Gene expression heatmap'

  xpr_htmp_res_tsv:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/xpr_htmp_res_tsv
    label: "Markers from gene expression heatmap"
    doc: |
      Gene markers used for gene
      expression heatmap

  gene_markers_tsv:
    type: File?
    outputSource: sc_rna_cluster/gene_markers_tsv
    label: "Gene markers per cluster for all resolutions"
    doc: |
      Gene markers per cluster for
      all resolutions
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Gene markers'
        Title: 'Gene markers per cluster for all resolutions'

  ucsc_cb_html_data:
    type: Directory?
    outputSource: sc_rna_cluster/ucsc_cb_html_data
    label: "UCSC Cell Browser data"
    doc: |
      Directory with UCSC Cell Browser
      data

  ucsc_cb_html_file:
    type: File?
    outputSource: sc_rna_cluster/ucsc_cb_html_file
    label: "UCSC Cell Browser"
    doc: |
      UCSC Cell Browser HTML index file
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  seurat_data_rds:
    type: File
    outputSource: sc_rna_cluster/seurat_data_rds
    label: "Processed Seurat data in RDS format"
    doc: |
      Processed Seurat data in RDS format

  seurat_data_scope:
    type: File?
    outputSource: sc_rna_cluster/seurat_data_scope
    label: "Processed Seurat data in SCope compatible loom format"
    doc: |
      Processed Seurat data in SCope compatible loom format

  pdf_plots:
    type: File
    outputSource: compress_pdf_plots/compressed_folder
    label: "Plots in PDF format"
    doc: |
      Compressed folder with plots
      in PDF format

  sc_rna_cluster_stdout_log:
    type: File
    outputSource: sc_rna_cluster/stdout_log
    label: "stdout log generated by sc_rna_cluster step"
    doc: |
      stdout log generated by sc_rna_cluster step

  sc_rna_cluster_stderr_log:
    type: File
    outputSource: sc_rna_cluster/stderr_log
    label: "stderr log generated by sc_rna_cluster step"
    doc: |
      stderr log generated by sc_rna_cluster step


steps:

  sc_rna_cluster:
    doc: |
      Clusters single-cell RNA-Seq datasets, identifies gene markers
    run: ../tools/sc-rna-cluster.cwl
    in:
      query_data_rds: query_data_rds
      dimensions: dimensions
      cluster_metric:
        default: euclidean
      cluster_algorithm:
        default: "louvain"
      resolution: resolution
      genes_of_interest:
        source: genes_of_interest
        valueFrom: $(split_features(self))
      identify_diff_genes: identify_diff_genes
      only_positive_diff_genes:
        default: true
      test_to_use: 
        default: wilcox
      minimum_logfc:
        default: 0.25
      minimum_pct:
        default: 0.1
      verbose:
        default: true
      export_ucsc_cb:
        default: true
      export_scope_data:
        default: true
      export_pdf_plots:
        default: true
      color_theme: color_theme
      parallel_memory_limit:
        default: 32
      vector_memory_limit:
        default: 96
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - umap_res_plot_png
    - slh_res_plot_png
    - umap_spl_idnt_res_plot_png
    - cmp_gr_clst_spl_idnt_res_plot_png
    - cmp_gr_idnt_spl_clst_res_plot_png
    - umap_spl_cnd_res_plot_png
    - cmp_gr_clst_spl_cnd_res_plot_png
    - cmp_gr_cnd_spl_clst_res_plot_png
    - umap_spl_ph_res_plot_png
    - cmp_gr_ph_spl_idnt_plot_png
    - cmp_gr_ph_spl_clst_res_plot_png
    - xpr_avg_res_plot_png
    - xpr_per_cell_plot_png
    - xpr_dnst_res_plot_png
    - xpr_htmp_res_plot_png
    - umap_res_plot_pdf
    - slh_res_plot_pdf
    - umap_spl_idnt_res_plot_pdf
    - cmp_gr_clst_spl_idnt_res_plot_pdf
    - cmp_gr_idnt_spl_clst_res_plot_pdf
    - umap_spl_cnd_res_plot_pdf
    - cmp_gr_clst_spl_cnd_res_plot_pdf
    - cmp_gr_cnd_spl_clst_res_plot_pdf
    - umap_spl_ph_res_plot_pdf
    - cmp_gr_ph_spl_idnt_plot_pdf
    - cmp_gr_ph_spl_clst_res_plot_pdf
    - xpr_avg_res_plot_pdf
    - xpr_per_cell_plot_pdf
    - xpr_per_cell_sgnl_plot_pdf
    - xpr_dnst_res_plot_pdf
    - xpr_htmp_res_plot_pdf
    - xpr_htmp_res_tsv
    - gene_markers_tsv
    - ucsc_cb_html_data
    - ucsc_cb_html_file
    - seurat_data_rds
    - seurat_data_scope
    - stdout_log
    - stderr_log

  pdf_plots:
    run: ../tools/files-to-folder.cwl
    in:
      input_files:
        source:
        - sc_rna_cluster/umap_res_plot_pdf
        - sc_rna_cluster/slh_res_plot_pdf
        - sc_rna_cluster/umap_spl_idnt_res_plot_pdf
        - sc_rna_cluster/cmp_gr_clst_spl_idnt_res_plot_pdf
        - sc_rna_cluster/cmp_gr_idnt_spl_clst_res_plot_pdf
        - sc_rna_cluster/umap_spl_cnd_res_plot_pdf
        - sc_rna_cluster/cmp_gr_clst_spl_cnd_res_plot_pdf
        - sc_rna_cluster/cmp_gr_cnd_spl_clst_res_plot_pdf
        - sc_rna_cluster/umap_spl_ph_res_plot_pdf
        - sc_rna_cluster/cmp_gr_ph_spl_idnt_plot_pdf
        - sc_rna_cluster/cmp_gr_ph_spl_clst_res_plot_pdf
        - sc_rna_cluster/xpr_avg_res_plot_pdf
        - sc_rna_cluster/xpr_per_cell_plot_pdf
        - sc_rna_cluster/xpr_per_cell_sgnl_plot_pdf
        - sc_rna_cluster/xpr_dnst_res_plot_pdf
        - sc_rna_cluster/xpr_htmp_res_plot_pdf
        valueFrom: $(self.flat().filter(n => n))
      folder_basename:
        default: "pdf_plots"
    out:
    - folder

  compress_pdf_plots:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: pdf_plots/folder
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-Cell RNA-Seq Cluster Analysis"
s:name: "Single-Cell RNA-Seq Cluster Analysis"
s:alternateName: "Clusters cells by similarity of gene expression data"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-rna-cluster.cwl
s:codeRepository: https://github.com/Barski-lab/workflows-datirium
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
  Single-Cell RNA-Seq Cluster Analysis

  Clusters cells by similarity of gene expression data from
  the outputs of “Single-Cell RNA-Seq Dimensionality Reduction
  Analysis” pipeline. The results of this workflow are primarily
  used in “Single-Cell Manual Cell Type Assignment” pipeline.