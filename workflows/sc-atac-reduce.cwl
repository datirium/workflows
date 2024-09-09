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

"sd:upstream":
  sc_tools_sample:
  - "sc-rna-cluster.cwl"
  - "sc-rna-reduce.cwl"
  - "sc-atac-filter.cwl"
  - "sc-multiome-filter.cwl"
  - "sc-rna-azimuth.cwl"


inputs:

  alias:
    type: string
    label: "Analysis name"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Single-cell Analysis with Filtered ATAC-Seq Datasets"
    doc: |
      Any analysis that includes single-cell
      multiome ATAC and RNA-Seq or just
      ATAC-Seq datasets filtered by QC metrics
      to include only high-quality cells.
    "sd:upstreamSource": "sc_tools_sample/seurat_data_rds"
    "sd:localLabel": true

  normalization_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "log-tfidf"
      - "tf-logidf"
      - "logtf-logidf"
      - "idf"
    label: "Normalization method"
    default: "log-tfidf"
    doc: |
      TF-IDF normalization method to correct
      for differences in cellular sequencing
      depth. "log-tfidf" - Stuart & Butler
      et al. 2019. "tf-logidf" - Cusanovich &
      Hill et al. 2018. "logtf-logidf" - Andrew
      Hill. "idf" - 10x Genomics. For more
      details refer to
      https://stuartlab.org/signac/reference/runtfidf
      Default: log-tfidf

  integration_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "signac"
      - "harmony"
      - "none"
    label: "Integration method"
    default: "signac"
    doc: |
      Integration method to match shared cell
      types and states across experimental
      batches, donors, conditions, or datasets.
      "signac" - use cross-dataset pairs of
      cells that are in a matched biological
      state ("anchors") to correct for technical
      differences. "harmony" - use Harmony
      algorithm described in Korsunsky, Millard,
      and Fan, Nat Methods, 2019, to iteratively
      correct LSI embeddings. "none" - do not
      run integration, merge datasets instead.
      Default: signac

  integrate_by:
    type:
    - "null"
    - string
    - type: enum
      symbols:
      - "dataset"
      - "condition"
    label: "Batch correction (harmony)"
    default: "dataset"
    doc: |
      When "harmony" is selected as "Integration
      method", batch effects are corrected based
      on the provided factors. Specifically,
      "dataset" is used to integrate out the
      influence of the cells' dataset of origin,
      while the factor "condition" is used to
      eliminate the influence of dataset grouping.
      Default: dataset

  dimensions:
    type: int?
    label: "Target dimensionality"
    default: 40
    doc: |
      Number of dimensions to be used in LSI,
      datasets integration, and UMAP projection.
      Accepted values range from 2 to 50.
      Default: 40

  datasets_metadata:
    type: File?
    label: "Datasets metadata (optional)"
    doc: |
      If the selected single-cell analysis
      includes multiple aggregated datasets,
      each of them can be assigned to a
      separate group by one or multiple
      categories. This can be achieved by
      providing a TSV/CSV file with
      "library_id" as the first column and
      any number of additional columns with
      unique names, representing the desired
      grouping categories. To obtain a proper
      template of this file, download
      "datasets_metadata.tsv" output from the
      "Files" tab of the selected "Single-cell
      Analysis with Filtered ATAC-Seq Datasets"
      and add extra columns as needed.

  barcodes_data:
    type: File?
    label: "Selected cell barcodes (optional)"
    doc: |
      A TSV/CSV file to optionally prefilter
      the single cell data by including only
      the cells with the selected barcodes.
      The provided file should include at
      least one column named "barcode", with
      one cell barcode per line. All other
      columns, except for "barcode", will be
      added to the single cell metadata loaded
      from "Single-cell Analysis with Filtered
      ATAC-Seq Datasets" and can be utilized in
      the current or future steps of analysis.

  minimum_var_peaks_perc:
    type: int?
    label: "Minimum percentile of highly variable peaks"
    default: 0
    doc: |
      Minimum percentile for identifying
      the top most common peaks as highly
      variable. For example, setting to 5
      will use the the top 95 percent most
      common among all cells peaks as highly
      variable. Selected peaks are then being
      used for datasets integration, scaling
      and dimensionality reduction.
      Default: 0 (use all available peaks)
    "sd:layout":
      advanced: true

  export_ucsc_cb:
    type: boolean?
    default: false
    label: "Show results in UCSC Cell Browser"
    doc: |
      Export results into UCSC Cell Browser
      Default: false
    "sd:layout":
      advanced: true

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

  qc_dim_corr_plot_png:
    type: File?
    outputSource: sc_atac_reduce/qc_dim_corr_plot_png
    label: "Correlation between QC metrics and LSI components"
    doc: |
      Correlation between QC metrics
      and LSI components
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Correlation between QC metrics and LSI components"

  umap_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_qc_mtrcs_plot_png
    label: "UMAP, QC metrics"
    doc: |
      UMAP, QC metrics
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "UMAP, QC metrics"

  umap_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_plot_png
    label: "UMAP, colored by dataset"
    doc: |
      UMAP, colored by dataset
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "UMAP, colored by dataset"

  umap_spl_idnt_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_spl_idnt_plot_png
    label: "UMAP, split by dataset"
    doc: |
      UMAP, split by dataset
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "UMAP, split by dataset"

  umap_spl_frgm_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_spl_frgm_plot_png
    label: "UMAP, colored by dataset, split by ATAC fragments in peaks per cell"
    doc: |
      UMAP, colored by dataset, split
      by ATAC fragments in peaks per cell.
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "UMAP, colored by dataset, split by ATAC fragments in peaks per cell"

  umap_spl_peak_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_spl_peak_plot_png
    label: "UMAP, colored by dataset, split by peaks per cell"
    doc: |
      UMAP, colored by dataset, split
      by peaks per cell
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "UMAP, colored by dataset, split by peaks per cell"

  umap_spl_tss_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_spl_tss_plot_png
    label: "UMAP, colored by dataset, split by TSS enrichment score"
    doc: |
      UMAP, colored by dataset, split
      by TSS enrichment score
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "UMAP, colored by dataset, split by TSS enrichment score"

  umap_spl_ncls_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_spl_ncls_plot_png
    label: "UMAP, colored by dataset, split by nucleosome signal"
    doc: |
      UMAP, colored by dataset, split
      by nucleosome signal
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "UMAP, colored by dataset, split by nucleosome signal"

  umap_spl_frip_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_spl_frip_plot_png
    label: "UMAP, colored by dataset, split by FRiP"
    doc: |
      UMAP, colored by dataset,
      split by FRiP
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "UMAP, colored by dataset, split by FRiP"

  umap_spl_blck_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_spl_blck_plot_png
    label: "UMAP, colored by dataset, split by blacklist fraction"
    doc: |
      UMAP, colored by dataset, split
      by blacklist fraction
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "UMAP, colored by dataset, split by blacklist fraction"

  umap_spl_cnd_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_spl_cnd_plot_png
    label: "UMAP, colored by dataset, split by grouping condition"
    doc: |
      UMAP, colored by dataset, split
      by grouping condition
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "UMAP, colored by dataset, split by grouping condition"

  umap_gr_cnd_spl_frgm_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_gr_cnd_spl_frgm_plot_png
    label: "UMAP, colored by grouping condition, split by ATAC fragments in peaks per cell"
    doc: |
      UMAP, colored by grouping condition,
      split by ATAC fragments in peaks per cell
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "UMAP, colored by grouping condition, split by ATAC fragments in peaks per cell"

  umap_gr_cnd_spl_peak_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_gr_cnd_spl_peak_plot_png
    label: "UMAP, colored by grouping condition, split by peaks per cell"
    doc: |
      UMAP, colored by grouping condition,
      split by peaks per cell
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "UMAP, colored by grouping condition, split by peaks per cell"

  umap_gr_cnd_spl_tss_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_gr_cnd_spl_tss_plot_png
    label: "UMAP, colored by grouping condition, split by TSS enrichment score"
    doc: |
      UMAP, colored by grouping condition,
      split by TSS enrichment score
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "UMAP, colored by grouping condition, split by TSS enrichment score"

  umap_gr_cnd_spl_ncls_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_gr_cnd_spl_ncls_plot_png
    label: "UMAP, colored by grouping condition, split by nucleosome signal"
    doc: |
      UMAP, colored by grouping condition,
      split by nucleosome signal
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "UMAP, colored by grouping condition, split by nucleosome signal"

  umap_gr_cnd_spl_frip_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_gr_cnd_spl_frip_plot_png
    label: "UMAP, colored by grouping condition, split by FRiP"
    doc: |
      UMAP, colored by grouping condition,
      split by FRiP
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "UMAP, colored by grouping condition, split by FRiP"

  umap_gr_cnd_spl_blck_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_gr_cnd_spl_blck_plot_png
    label: "UMAP, colored by grouping condition, split by blacklist fraction"
    doc: |
      UMAP, colored by grouping condition,
      split by blacklist fraction
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "UMAP, colored by grouping condition, split by blacklist fraction"

  ucsc_cb_html_data:
    type: Directory?
    outputSource: sc_atac_reduce/ucsc_cb_html_data
    label: "UCSC Cell Browser (data)"
    doc: |
      UCSC Cell Browser html data.

  ucsc_cb_html_file:
    type: File?
    outputSource: sc_atac_reduce/ucsc_cb_html_file
    label: "UCSC Cell Browser"
    doc: |
      UCSC Cell Browser html index.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  seurat_data_rds:
    type: File
    outputSource: sc_atac_reduce/seurat_data_rds
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
    outputSource: sc_atac_reduce/sc_report_html_file
    label: "Analysis log"
    doc: |
      Tehcnical report.
      HTML format.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  sc_atac_reduce_stdout_log:
    type: File
    outputSource: sc_atac_reduce/stdout_log
    label: "Output log"
    doc: |
      Stdout log from the sc_atac_reduce step.

  sc_atac_reduce_stderr_log:
    type: File
    outputSource: sc_atac_reduce/stderr_log
    label: "Error log"
    doc: |
      Stderr log from the sc_atac_reduce step.


steps:

  sc_atac_reduce:
    doc: |
      Integrates multiple single-cell ATAC-Seq datasets,
      reduces dimensionality using LSI
    run: ../tools/sc-atac-reduce.cwl
    in:
      query_data_rds: query_data_rds
      barcodes_data: barcodes_data
      datasets_metadata: datasets_metadata
      normalization_method: normalization_method
      integration_method: integration_method
      integrate_by:
        source: integrate_by
        valueFrom: |
          ${
            if (self == "none") {
              return null;
            } else if (self == "dataset") {
              return "new.ident";
            } else if (self == "condition") {
              return "condition";
            } else {
              return split_features(self);
            }
          }
      minimum_var_peaks_perc: minimum_var_peaks_perc
      dimensions: dimensions
      verbose:
        default: true
      export_ucsc_cb: export_ucsc_cb
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
    - qc_dim_corr_plot_png
    - umap_qc_mtrcs_plot_png
    - umap_plot_png
    - umap_spl_idnt_plot_png
    - umap_spl_cnd_plot_png
    - umap_spl_frgm_plot_png
    - umap_spl_peak_plot_png
    - umap_spl_tss_plot_png
    - umap_spl_ncls_plot_png
    - umap_spl_frip_plot_png
    - umap_spl_blck_plot_png
    - umap_gr_cnd_spl_frgm_plot_png
    - umap_gr_cnd_spl_peak_plot_png
    - umap_gr_cnd_spl_tss_plot_png
    - umap_gr_cnd_spl_ncls_plot_png
    - umap_gr_cnd_spl_frip_plot_png
    - umap_gr_cnd_spl_blck_plot_png
    - qc_dim_corr_plot_pdf
    - umap_qc_mtrcs_plot_pdf
    - umap_plot_pdf
    - umap_spl_idnt_plot_pdf
    - umap_spl_cnd_plot_pdf
    - umap_spl_frgm_plot_pdf
    - umap_spl_peak_plot_pdf
    - umap_spl_tss_plot_pdf
    - umap_spl_ncls_plot_pdf
    - umap_spl_frip_plot_pdf
    - umap_spl_blck_plot_pdf
    - umap_gr_cnd_spl_frgm_plot_pdf
    - umap_gr_cnd_spl_peak_plot_pdf
    - umap_gr_cnd_spl_tss_plot_pdf
    - umap_gr_cnd_spl_ncls_plot_pdf
    - umap_gr_cnd_spl_frip_plot_pdf
    - umap_gr_cnd_spl_blck_plot_pdf
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
        - sc_atac_reduce/qc_dim_corr_plot_pdf
        - sc_atac_reduce/umap_qc_mtrcs_plot_pdf
        - sc_atac_reduce/umap_plot_pdf
        - sc_atac_reduce/umap_spl_idnt_plot_pdf
        - sc_atac_reduce/umap_spl_cnd_plot_pdf
        - sc_atac_reduce/umap_spl_frgm_plot_pdf
        - sc_atac_reduce/umap_spl_peak_plot_pdf
        - sc_atac_reduce/umap_spl_tss_plot_pdf
        - sc_atac_reduce/umap_spl_ncls_plot_pdf
        - sc_atac_reduce/umap_spl_frip_plot_pdf
        - sc_atac_reduce/umap_spl_blck_plot_pdf
        - sc_atac_reduce/umap_gr_cnd_spl_frgm_plot_pdf
        - sc_atac_reduce/umap_gr_cnd_spl_peak_plot_pdf
        - sc_atac_reduce/umap_gr_cnd_spl_tss_plot_pdf
        - sc_atac_reduce/umap_gr_cnd_spl_ncls_plot_pdf
        - sc_atac_reduce/umap_gr_cnd_spl_frip_plot_pdf
        - sc_atac_reduce/umap_gr_cnd_spl_blck_plot_pdf
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

label: "Single-Cell ATAC-Seq Dimensionality Reduction Analysis"
s:name: "Single-Cell ATAC-Seq Dimensionality Reduction Analysis"
s:alternateName: "Removes noise and confounding sources of variation by reducing dimensionality of chromatin accessibility data"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-atac-reduce.cwl
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
  Single-Cell ATAC-Seq Dimensionality Reduction Analysis

  Removes noise and confounding sources of variation by reducing
  dimensionality of chromatin accessibility data from the outputs
  of “Single-Cell Multiome ATAC and RNA-Seq Filtering Analysis”
  pipelines. The results of this workflow are primarily used in
  “Single-Cell ATAC-Seq Cluster Analysis” or “Single-Cell WNN
  Cluster Analysis” pipelines.