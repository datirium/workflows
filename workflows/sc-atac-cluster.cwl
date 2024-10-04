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
    - var parse_range = function(line) {
          if (line.includes("-")) {
              const parts = line.split("-");
              const start = parseFloat(parts[0].trim());
              let end, step;
              if (parts[1].includes(":")) {
                  [end, step] = parts[1].split(":").map(Number);
              } else {
                  end = parseFloat(parts[1].trim());
                  step = 0.1;
              }
              const result = [];
              for (let i = start; i <= end; i = parseFloat((i + step).toFixed(10))) {
                  result.push(parseFloat(i.toFixed(10)));
              }
              return result;
          } else {
              return [parseFloat(line)];
          }
      };


"sd:upstream":
  sc_tools_sample:
  - "sc-atac-cluster.cwl"
  - "sc-rna-cluster.cwl"
  - "sc-rna-reduce.cwl"
  - "sc-atac-reduce.cwl"
  - "sc-rna-azimuth.cwl"
  sc_atac_sample:
  - "cellranger-arc-count.cwl"
  - "cellranger-arc-aggr.cwl"
  - "cellranger-atac-count.cwl"
  - "cellranger-atac-aggr.cwl"


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
      "Single-Cell ATAC-Seq Dimensionality
      Reduction Analysis" at any of the
      processing stages.
    "sd:upstreamSource": "sc_tools_sample/seurat_data_rds"
    "sd:localLabel": true

  atac_fragments_file:
    type: File?
    secondaryFiles:
    - .tbi
    label: "Cell Ranger ATAC or RNA+ATAC Sample (optional)"
    doc: |
      Any "Cell Ranger ATAC or RNA+ATAC Sample"
      for generating ATAC fragments coverage
      plots over the genes of interest. This
      sample can be obtained from one of the
      following pipelines: "Cell Ranger Count
      (RNA+ATAC)", "Cell Ranger Aggregate
      (RNA+ATAC)", "Cell Ranger Count (ATAC)",
      or "Cell Ranger Aggregate (ATAC)".
    "sd:upstreamSource": "sc_atac_sample/atac_fragments_file"
    "sd:localLabel": true

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
    type: string?
    default: "0.3"
    label: "Clustering resolution"
    doc: |
      Resolution to define the "granularity"
      of the clustered data. Larger values
      lead to a bigger number of clusters.
      Optimal resolution often increases
      with the number of cells. To run the
      analysis with multiple resolutions,
      provide a range in a form of
      start-end:step. Step parameter is
      optional and equal to 0.1 by default.
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
      of interest to generate ATAC fragments
      coverage plots. Ignored if "Cell Ranger
      ATAC or RNA+ATAC (optional)" input is
      not provided.
      Default: None

  export_html_report:
    type: boolean?
    default: true
    label: "Show HTML report"
    doc: |
      Export tehcnical report in HTML format.
      Default: true
    "sd:layout":
      advanced: true

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
    default: "4"
    label: "Cores/CPUs"
    doc: |
      Parallelization parameter to define the
      number of cores/CPUs that can be utilized
      simultaneously.
      Default: 4
    "sd:layout":
      advanced: true


outputs:

  cell_cnts_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_cluster/cell_cnts_gr_clst_res_plot_png
    label: "Number of cells per cluster (all cells)"
    doc: |
      Number of cells per cluster.
      All cells; all resolutions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Number of cells per cluster (all cells)"

  qc_mtrcs_dnst_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_cluster/qc_mtrcs_dnst_gr_clst_res_plot_png
    label: "Distribution of QC metrics per cell colored by cluster (all cells)"
    doc: |
      Distribution of QC metrics per cell
      colored by cluster.
      All cells; all resolutions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Distribution of QC metrics per cell colored by cluster (all cells)"

  tss_frgm_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_cluster/tss_frgm_spl_clst_res_plot_png
    label: "TSS enrichment score vs ATAC fragments in peaks per cell (split by cluster, all cells)"
    doc: |
      TSS enrichment score vs ATAC
      fragments in peaks per cell.
      Split by cluster; all cells;
      all resolutions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "TSS enrichment score vs ATAC fragments in peaks per cell (split by cluster, all cells)"

  atacdbl_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_cluster/atacdbl_gr_clst_res_plot_png
    label: "Percentage of ATAC doublets per cluster (all cells)"
    doc: |
      Percentage of ATAC doublets per cluster.
      All cells; all resolutions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Percentage of ATAC doublets per cluster (all cells)"

  umap_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_cluster/umap_gr_clst_res_plot_png
    label: "UMAP colored by cluster (all cells)"
    doc: |
      UMAP colored by cluster.
      All cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by cluster"
        Caption: "UMAP colored by cluster (all cells)"

  slh_gr_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_cluster/slh_gr_clst_res_plot_png
    label: "Silhouette scores (all cells)"
    doc: |
      Silhouette scores.
      All cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by cluster"
        Caption: "Silhouette scores (all cells)"

  umap_gr_clst_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_cluster/umap_gr_clst_spl_idnt_res_plot_png
    label: "UMAP colored by cluster (split by dataset, downsampled)"
    doc: |
      UMAP colored by cluster.
      Split by dataset; downsampled
      to the smallest dataset.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by dataset"
        Caption: "UMAP colored by cluster (split by dataset, downsampled)"

  cmp_gr_clst_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_cluster/cmp_gr_clst_spl_idnt_res_plot_png
    label: "Composition plot colored by cluster (split by dataset, downsampled)"
    doc: |
      Composition plot colored by cluster.
      Split by dataset; downsampled
      to the smallest dataset.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by dataset"
        Caption: "Composition plot colored by cluster (split by dataset, downsampled)"

  umap_gr_clst_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_cluster/umap_gr_clst_spl_cnd_res_plot_png
    label: "UMAP colored by cluster (split by grouping condition, downsampled)"
    doc: |
      UMAP colored by cluster.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by group"
        Caption: "UMAP colored by cluster (split by grouping condition, downsampled)"

  cmp_gr_clst_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_cluster/cmp_gr_clst_spl_cnd_res_plot_png
    label: "Composition plot colored by cluster (split by grouping condition, downsampled)"
    doc: |
      Composition plot colored by cluster.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by group"
        Caption: "Composition plot colored by cluster (split by grouping condition, downsampled)"

  cvrg_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_cluster/cvrg_res_plot_png
    label: "ATAC fragment coverage (per gene)"
    doc: |
      ATAC fragment coverage.
      All genes of interest.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Genes of interest (coverage plot)"
        Caption: "ATAC fragment coverage (per gene)"

  peak_markers_tsv:
    type: File?
    outputSource: sc_atac_cluster/peak_markers_tsv
    label: "Peak markers"
    doc: |
      Peak markers.
      TSV format.
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Peak markers table"
        Title: "Peak markers"

  ucsc_cb_html_data:
    type: Directory?
    outputSource: sc_atac_cluster/ucsc_cb_html_data
    label: "UCSC Cell Browser (data)"
    doc: |
      UCSC Cell Browser html data.

  ucsc_cb_html_file:
    type: File?
    outputSource: sc_atac_cluster/ucsc_cb_html_file
    label: "UCSC Cell Browser"
    doc: |
      UCSC Cell Browser html index.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  seurat_data_rds:
    type: File
    outputSource: sc_atac_cluster/seurat_data_rds
    label: "Seurat object in RDS format"
    doc: |
      Seurat object.
      RDS format.

  pdf_plots:
    type: File
    outputSource: compress_pdf_plots/compressed_folder
    label: "Compressed folder with all PDF plots"
    doc: |
      Compressed folder with all PDF plots.

  sc_report_html_file:
    type: File?
    outputSource: sc_atac_cluster/sc_report_html_file
    label: "Analysis log"
    doc: |
      Tehcnical report.
      HTML format.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  sc_atac_cluster_stdout_log:
    type: File
    outputSource: sc_atac_cluster/stdout_log
    label: "Output log"
    doc: |
      Stdout log from the sc_atac_cluster step.

  sc_atac_cluster_stderr_log:
    type: File
    outputSource: sc_atac_cluster/stderr_log
    label: "Error log"
    doc: |
      Stderr log from the sc_atac_cluster step.


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
      resolution:
        source: resolution
        valueFrom: $(parse_range(self))
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
      export_pdf_plots:
        default: true
      color_theme: color_theme
      parallel_memory_limit:
        default: 32
      vector_memory_limit:
        default: 128
      export_html_report: export_html_report
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - cell_cnts_gr_clst_res_plot_png
    - tss_frgm_spl_clst_res_plot_png
    - atacdbl_gr_clst_res_plot_png
    - qc_mtrcs_dnst_gr_clst_res_plot_png
    - umap_gr_clst_res_plot_png
    - slh_gr_clst_res_plot_png
    - umap_gr_clst_spl_idnt_res_plot_png
    - cmp_gr_clst_spl_idnt_res_plot_png
    - umap_gr_clst_spl_cnd_res_plot_png
    - cmp_gr_clst_spl_cnd_res_plot_png
    - cvrg_res_plot_png
    - all_plots_pdf
    - peak_markers_tsv
    - ucsc_cb_html_data
    - ucsc_cb_html_file
    - seurat_data_rds
    - sc_report_html_file
    - stdout_log
    - stderr_log

  folder_pdf_plots:
    run: ../tools/files-to-folder.cwl
    in:
      input_files:
        source:
        - sc_atac_cluster/all_plots_pdf
        valueFrom: $(self.flat().filter(n => n))
      folder_basename:
        default: "pdf_plots"
    out:
    - folder

  compress_pdf_plots:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: folder_pdf_plots/folder
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-Cell ATAC-Seq Cluster Analysis"
s:name: "Single-Cell ATAC-Seq Cluster Analysis"
s:alternateName: "Single-Cell ATAC-Seq Cluster Analysis"

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
  Single-Cell ATAC-Seq Cluster Analysis

  Clusters cells by similarity of chromatin accessibility
  data from the outputs of the “Single-Cell ATAC-Seq
  Dimensionality Reduction Analysis” pipeline. The results
  of this workflow are used in the “Single-Cell Manual Cell
  Type Assignment”, “Single-Cell ATAC-Seq Differential
  Accessibility Analysis”, and “Single-Cell ATAC-Seq Genome
  Coverage” pipelines.