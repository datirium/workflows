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
  - "sc-atac-cluster.cwl"
  - "sc-rna-cluster.cwl"
  - "sc-rna-reduce.cwl"
  - "sc-atac-reduce.cwl"
  sc_arc_sample:
  - "cellranger-arc-count.cwl"
  - "cellranger-arc-aggr.cwl"
  - "https://github.com/datirium/workflows/workflows/cellranger-arc-count.cwl"
  - "https://github.com/datirium/workflows/workflows/cellranger-arc-aggr.cwl"


inputs:

  alias:
    type: string
    label: "Analysis name"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Single-cell Analysis with LSI Transformed ATAC-Seq Datasets"
    doc: |
      Analysis that includes single-cell
      multiome RNA and ATAC-Seq or just
      ATAC-Seq datasets run through
      "Single-cell ATAC-Seq Dimensionality
      Reduction Analysis" at any of the
      processing stages.
    'sd:upstreamSource': "sc_tools_sample/seurat_data_rds"
    'sd:localLabel': true

  atac_fragments_file:
    type: File?
    secondaryFiles:
    - .tbi
    label: "Cell Ranger ARC Sample (optional)"
    doc: |
      "Cell Ranger ARC Sample" for generating
      fragments coverage plots over the genes
      of interest.
    'sd:upstreamSource': "sc_arc_sample/atac_fragments_file"
    'sd:localLabel': true

  dimensions:
    type: int?
    default: 40
    label: "Target dimensionality"
    doc: |
      Number of LSI components to be used
      in constructing nearest-neighbor graph
      as part of the clustering algorithm.
      Accepted values range from 2 to 50.
      First dimension is always excluded
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
      with the number of cells.
      Default: 0.3

  identify_diff_peaks:
    type: boolean?
    default: false
    label: "Find peak markers"
    doc: |
      Identify differentially accessible
      peaks in each cluster compared to
      all other cells. Include only peaks
      that are present in at least 5% of
      the cells coming from either current
      cluster or from all other clusters
      together. Exclude cells with
      log2FoldChange values less than 0.25.
      Use logistic regression framework to
      calculate P-values. Keep only genes
      with P-values lower than 0.01. Adjust
      P-values for multiple comparisons
      using Bonferroni correction.
      Default: false

  genes_of_interest:
    type: string?
    default: null
    label: "Genes of interest"
    doc: |
      Comma or space separated list of genes
      of interest to generate fragments coverage
      plots. Ignored if "Cell Ranger ARC Sample"
      input is not provided.
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
    default: "1"
    label: "Cores/CPUs number"
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
    outputSource: sc_atac_cluster/umap_res_plot_png
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
    outputSource: sc_atac_cluster/slh_res_plot_png
    label: "Silhouette scores"
    doc: |
      Silhouette scores
    'sd:visualPlugins':
    - image:
        tab: 'Per cluster'
        Caption: 'Silhouette scores'

  umap_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_cluster/umap_spl_idnt_res_plot_png
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
    outputSource: sc_atac_cluster/cmp_gr_clst_spl_idnt_res_plot_png
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
    outputSource: sc_atac_cluster/cmp_gr_idnt_spl_clst_res_plot_png
    label: "Composition plot, colored by dataset, split by cluster, downsampled"
    doc: |
      Composition plot, colored by
      dataset, split by cluster,
      downsampled
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Composition plot, colored by dataset, split by cluster, downsampled'

  umap_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_cluster/umap_spl_cnd_res_plot_png
    label: "UMAP, colored by cluster, split by grouping condition"
    doc: |
      UMAP, colored by cluster, split
      by grouping condition
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'UMAP, colored by cluster, split by grouping condition'

  cmp_gr_clst_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_cluster/cmp_gr_clst_spl_cnd_res_plot_png
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
    outputSource: sc_atac_cluster/cmp_gr_cnd_spl_clst_res_plot_png
    label: "Composition plot, colored by grouping condition, split by cluster, downsampled"
    doc: |
      Composition plot, colored by
      grouping condition, split by
      cluster, downsampled
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Composition plot, colored by grouping condition, split by cluster, downsampled'

  cvrg_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_cluster/cvrg_res_plot_png
    label: "Fragments coverage"
    doc: |
      Fragments coverage
    'sd:visualPlugins':
    - image:
        tab: 'Genome coverage'
        Caption: 'Fragments coverage'

  peak_markers_tsv:
    type: File?
    outputSource: sc_atac_cluster/peak_markers_tsv
    label: "Peak markers per cluster for all resolutions"
    doc: |
      Peak markers per cluster for all resolutions
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Peak markers'
        Title: 'Peak markers per cluster for all resolutions'

  ucsc_cb_html_data:
    type: Directory?
    outputSource: sc_atac_cluster/ucsc_cb_html_data
    label: "UCSC Cell Browser data"
    doc: |
      Directory with UCSC Cell Browser
      data

  ucsc_cb_html_file:
    type: File?
    outputSource: sc_atac_cluster/ucsc_cb_html_file
    label: "UCSC Cell Browser"
    doc: |
      UCSC Cell Browser HTML index file
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  seurat_data_rds:
    type: File
    outputSource: sc_atac_cluster/seurat_data_rds
    label: "Processed Seurat data in RDS format"
    doc: |
      Processed Seurat data in RDS format

  sc_atac_cluster_stdout_log:
    type: File
    outputSource: sc_atac_cluster/stdout_log
    label: "stdout log generated by sc_atac_cluster step"
    doc: |
      stdout log generated by sc_atac_cluster step

  sc_atac_cluster_stderr_log:
    type: File
    outputSource: sc_atac_cluster/stderr_log
    label: "stderr log generated by sc_atac_cluster step"
    doc: |
      stderr log generated by sc_atac_cluster step


steps:

  sc_atac_cluster:
    doc: |
      Clusters single-cell ATAC-Seq datasets, identifies differentially
      accessible peaks
    run: ../tools/sc-atac-cluster.cwl
    in:
      query_data_rds: query_data_rds
      dimensions: dimensions
      cluster_metric:
        default: euclidean
      cluster_algorithm:
        default: "slm"
      resolution: resolution
      atac_fragments_file: atac_fragments_file
      genes_of_interest:
        source: genes_of_interest
        valueFrom: $(split_features(self))
      identify_diff_peaks: identify_diff_peaks
      minimum_logfc:
        default: 0.25
      minimum_pct:
        default: 0.05
      test_to_use: 
        default: LR
      verbose:
        default: true
      export_ucsc_cb:
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
    - cvrg_res_plot_png
    - peak_markers_tsv
    - ucsc_cb_html_data
    - ucsc_cb_html_file
    - seurat_data_rds
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell ATAC-Seq Cluster Analysis"
s:name: "Single-cell ATAC-Seq Cluster Analysis"
s:alternateName: "Clusters single-cell ATAC-Seq datasets, identifies differentially accessible peaks"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-atac-cluster.cwl
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
  Single-cell ATAC-Seq Cluster Analysis

  Clusters single-cell ATAC-Seq datasets, identifies
  differentially accessible peaks.