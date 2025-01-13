cwlVersion: v1.1
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
  - "sc-atac-cluster.cwl"
  - "sc-atac-reduce.cwl"
  - "sc-rna-filter.cwl"
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
    label: "Single-cell Analysis with Filtered RNA-Seq Datasets"
    doc: |
      Any analysis that includes single-cell
      multiome ATAC and RNA-Seq or just
      RNA-Seq datasets filtered by QC metrics
      to include only high-quality cells.
    "sd:upstreamSource": "sc_tools_sample/seurat_data_rds"
    "sd:localLabel": true

  normalization_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "sct"
      - "sctglm"
      - "log"
    label: "Normalization method"
    default: "sctglm"
    doc: |
      Normalization and variance stabilization
      method to remove technical variability
      between the cells. "sct" - use sctransform
      package described in Hafemeister and Satija,
      Genome Biology 2019. "sctglm" - use updated
      sctransform package described in Choudhary
      and Satija, Genome Biology, 2022. "log" -
      use a combination of NormalizeData and
      ScaleData functions described in Stuart and
      Butler, Cell 2019.
      Default: sctglm

  integration_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "seurat"
      - "harmony"
      - "none"
    label: "Integration method"
    default: "seurat"
    doc: |
      Integration method to match shared cell
      types and states across experimental
      batches, donors, conditions, or datasets.
      "seurat" - use cross-dataset pairs of
      cells that are in a matched biological
      state ("anchors") to correct for technical
      differences. "harmony" - use Harmony
      algorithm described in Korsunsky, Millard,
      and Fan, Nat Methods, 2019, to iteratively
      correct PCA embeddings. "none" - do not
      run integration, merge datasets instead.
      Default: seurat

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
      Number of principal components to be used
      in PCA and UMAP projection. Accepted values
      range from 1 to 50. Set to 0 to use
      auto-estimated dimensionality.
      Default: 40

  cell_cycle_data:
    type:
    - "null"
    - type: enum
      symbols:
      - "human"
      - "mouse"
      - "none"
    label: "Cell cycle gene set"
    default: "none"
    doc: |
      Assign cell cycle score and
      phase based on the gene set
      for the selected organism.
      When selected "none", skip
      cell cycle score assignment.
      Default: "none"

  regress_cellcycle:
    type:
    - "null"
    - type: enum
      symbols:
      - "completely"
      - "partially"
      - "do not remove"
    label: "Remove cell cycle"
    default: "do not remove"
    doc: |
      Remove the influence cell cycle
      phase on the dimensionality
      reduction results. When selected
      "completely", regress all signals
      associated with the cell cycle phase.
      For "partially" - regress only the
      differences in cell cycle phase
      among proliferating cells, signals
      separating non-cycling and cycling
      cells will be maintained. When
      selected "do not remove" - do not
      regress signals associated with the
      cell cycle phase. Ignored if cell
      cycle gene set is not provided.
      Default: "do not remove"

  regress_genes:
    type: string?
    label: "Regress genes"
    default: null
    doc: |
      Regex pattern to identify genes which
      expression should be regressed as a
      confounding source of variation.
      Default: None

  regress_mito_perc:
    type: boolean?
    label: "Regress mitochondrial percentage"
    default: false
    doc: |
      Regress the percentage of RNA reads
      mapped to mitochondrial genes as a
      confounding source of variation.
      Default: false

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
      Analysis with Filtered RNA-Seq Datasets"
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
      RNA-Seq Datasets" and can be utilized in
      the current or future steps of analysis.

  custom_cell_cycle_data:
    type: File?
    label: "Custom cell cycle gene set (optional)"
    doc: |
      A TSV/CSV file with the gene list
      for cell cycle score assignment.
      The file should have two columns
      named 'phase' and 'gene_id'. If
      this input is provided, the "Cell
      cycle gene set" will be ignored.

  highly_var_genes_count:
    type: int?
    label: "Number of highly variable genes"
    default: 3000
    doc: |
      The number of highly variable genes
      to be used in gene expression scaling,
      datasets integration, and dimensionality
      reduction.
      Default: 3000
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

  export_loupe_data:
    type: boolean?
    default: false
    label: "Save raw counts to Loupe file. I confirm that data is generated by 10x technology and accept the EULA available at https://10xgen.com/EULA"
    doc: |
      Save raw counts from the RNA assay to Loupe file. By
      enabling this feature you accept the End-User License
      Agreement available at https://10xgen.com/EULA.
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

  elbow_plot_png:
    type: File?
    outputSource: sc_rna_reduce/elbow_plot_png
    label: "Elbow plot"
    doc: |
      Elbow plot to evaluate the number of
      principal components that capture the
      majority of the variation in the data.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Elbow plot"

  qc_dim_corr_plot_png:
    type: File?
    outputSource: sc_rna_reduce/qc_dim_corr_plot_png
    label: "Correlation between QC metrics and principal components"
    doc: |
      Correlation between QC metrics and
      principal components
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Correlation between QC metrics and principal components"

  umap_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_qc_mtrcs_plot_png
    label: "UMAP, QC metrics"
    doc: |
      UMAP, QC metrics
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "UMAP, QC metrics"

  ccpca_plot_png:
    type: File?
    outputSource: sc_rna_reduce/ccpca_plot_png
    label: "PCA, colored by cell cycle phase"
    doc: |
      PCA, colored by cell cycle phase
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "PCA, colored by cell cycle phase"

  umap_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_plot_png
    label: "UMAP, colored by dataset"
    doc: |
      UMAP, colored by dataset
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "UMAP, colored by dataset"

  umap_spl_idnt_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_idnt_plot_png
    label: "UMAP, split by dataset"
    doc: |
      UMAP, split by dataset
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "UMAP, split by dataset"

  umap_spl_umi_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_umi_plot_png
    label: "UMAP, colored by dataset, split by RNA reads per cell"
    doc: |
      UMAP, colored by dataset, split by
      RNA reads per cell
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "UMAP, colored by dataset, split by RNA reads per cell"

  umap_spl_gene_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_gene_plot_png
    label: "UMAP, colored by dataset, split by genes per cell"
    doc: |
      UMAP, colored by dataset, split by
      genes per cell
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "UMAP, colored by dataset, split by genes per cell"

  umap_spl_mito_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_mito_plot_png
    label: "UMAP, colored by dataset, split by mitochondrial percentage"
    doc: |
      UMAP, colored by dataset, split by
      mitochondrial percentage
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "UMAP, colored by dataset, split by mitochondrial percentage"

  umap_spl_ph_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_ph_plot_png
    label: "UMAP, colored by dataset, split by cell cycle phase"
    doc: |
      UMAP, colored by dataset, split by
      cell cycle phase
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "UMAP, colored by dataset, split by cell cycle phase"

  ccpca_spl_idnt_plot_png:
    type: File?
    outputSource: sc_rna_reduce/ccpca_spl_idnt_plot_png
    label: "PCA, colored by cell cycle phase, split by dataset"
    doc: |
      PCA, colored by cell cycle phase,
      split by dataset
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "PCA, colored by cell cycle phase, split by dataset"

  umap_spl_cnd_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_cnd_plot_png
    label: "UMAP, colored by dataset, split by grouping condition"
    doc: |
      UMAP, colored by dataset, split by
      grouping condition
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "UMAP, colored by dataset, split by grouping condition"

  umap_gr_cnd_spl_umi_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_umi_plot_png
    label: "UMAP, colored by grouping condition, split by RNA reads per cell"
    doc: |
      UMAP, colored by grouping condition,
      split by RNA reads per cell
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "UMAP, colored by grouping condition, split by RNA reads per cell"

  umap_gr_cnd_spl_gene_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_gene_plot_png
    label: "UMAP, colored by grouping condition, split by genes per cell"
    doc: |
      UMAP, colored by grouping condition,
      split by genes per cell
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "UMAP, colored by grouping condition, split by genes per cell"

  umap_gr_cnd_spl_mito_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_mito_plot_png
    label: "UMAP, colored by grouping condition, split by mitochondrial percentage"
    doc: |
      UMAP, colored by grouping condition,
      split by mitochondrial percentage
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "UMAP, colored by grouping condition, split by mitochondrial percentage"

  umap_gr_cnd_spl_ph_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_ph_plot_png
    label: "UMAP, colored by grouping condition, split by cell cycle phase"
    doc: |
      UMAP, colored by grouping condition,
      split by cell cycle phase
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "UMAP, colored by grouping condition, split by cell cycle phase"

  ccpca_spl_cnd_plot_png:
    type: File?
    outputSource: sc_rna_reduce/ccpca_spl_cnd_plot_png
    label: "PCA, colored by cell cycle phase, split by grouping condition"
    doc: |
      PCA, colored by cell cycle phase,
      split by grouping condition
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "PCA, colored by cell cycle phase, split by grouping condition"

  ucsc_cb_html_data:
    type: Directory?
    outputSource: sc_rna_reduce/ucsc_cb_html_data
    label: "UCSC Cell Browser (data)"
    doc: |
      UCSC Cell Browser html data.

  ucsc_cb_html_file:
    type: File?
    outputSource: sc_rna_reduce/ucsc_cb_html_file
    label: "UCSC Cell Browser"
    doc: |
      UCSC Cell Browser html index.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  seurat_data_rds:
    type: File
    outputSource: sc_rna_reduce/seurat_data_rds
    label: "Seurat object in RDS format"
    doc: |
      Seurat object.
      RDS format.

  seurat_data_cloupe:
    type: File?
    outputSource: sc_rna_reduce/seurat_data_cloupe
    label: "Seurat object in Loupe format"
    doc: |
      Seurat object.
      Loupe format.

  pdf_plots:
    type: File
    outputSource: compress_pdf_plots/compressed_folder
    label: "Compressed folder with all PDF plots"
    doc: |
      Compressed folder with all PDF plots.

  sc_report_html_file:
    type: File?
    outputSource: sc_rna_reduce/sc_report_html_file
    label: "Analysis log"
    doc: |
      Tehcnical report.
      HTML format.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  sc_rna_reduce_stdout_log:
    type: File
    outputSource: sc_rna_reduce/stdout_log
    label: "Output log"
    doc: |
      Stdout log from the sc_rna_reduce step.

  sc_rna_reduce_stderr_log:
    type: File
    outputSource: sc_rna_reduce/stderr_log
    label: "Error log"
    doc: |
      Stderr log from the sc_rna_reduce step.


steps:

  sc_rna_reduce:
    doc: |
      Integrates multiple single-cell RNA-Seq datasets,
      reduces dimensionality using PCA
    run: ../tools/sc-rna-reduce.cwl
    in:
      query_data_rds: query_data_rds
      barcodes_data: barcodes_data
      cell_cycle_data:
        source: [cell_cycle_data, custom_cell_cycle_data]
        valueFrom: |
          ${
            if (self[1] != null && self[1].class == "File"){
              return self[1];
            } else if (self[0].includes("human")) {
              return "hg38";
            } else if (self[0].includes("mouse")) {
              return "mm10";
            } else {
              return null;
            }
          }
      regress_ccycle_full:
        source: regress_cellcycle
        valueFrom: $(self.includes("completely")?true:null)
      regress_ccycle_diff:
        source: regress_cellcycle
        valueFrom: $(self.includes("partially")?true:null)
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
      highly_var_genes_count: highly_var_genes_count
      regress_mito_perc: regress_mito_perc
      regress_genes:
        source: regress_genes
        valueFrom: $(self==""?null:self)            # safety measure
      dimensions: dimensions
      verbose:
        default: true
      export_ucsc_cb: export_ucsc_cb
      low_memory:
        default: true
      export_loupe_data: export_loupe_data
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
    - elbow_plot_png
    - qc_dim_corr_plot_png
    - umap_qc_mtrcs_plot_png
    - umap_plot_png
    - umap_spl_ph_plot_png
    - ccpca_plot_png
    - umap_spl_mito_plot_png
    - umap_spl_umi_plot_png
    - umap_spl_gene_plot_png
    - umap_spl_idnt_plot_png
    - ccpca_spl_idnt_plot_png
    - umap_spl_cnd_plot_png
    - umap_gr_cnd_spl_ph_plot_png
    - ccpca_spl_cnd_plot_png
    - umap_gr_cnd_spl_mito_plot_png
    - umap_gr_cnd_spl_umi_plot_png
    - umap_gr_cnd_spl_gene_plot_png
    - elbow_plot_pdf
    - qc_dim_corr_plot_pdf
    - umap_qc_mtrcs_plot_pdf
    - umap_plot_pdf
    - umap_spl_ph_plot_pdf
    - ccpca_plot_pdf
    - umap_spl_mito_plot_pdf
    - umap_spl_umi_plot_pdf
    - umap_spl_gene_plot_pdf
    - umap_spl_idnt_plot_pdf
    - ccpca_spl_idnt_plot_pdf
    - umap_spl_cnd_plot_pdf
    - umap_gr_cnd_spl_ph_plot_pdf
    - ccpca_spl_cnd_plot_pdf
    - umap_gr_cnd_spl_mito_plot_pdf
    - umap_gr_cnd_spl_umi_plot_pdf
    - umap_gr_cnd_spl_gene_plot_pdf
    - ucsc_cb_html_data
    - ucsc_cb_html_file
    - seurat_data_rds
    - seurat_data_cloupe
    - sc_report_html_file
    - stdout_log
    - stderr_log

  folder_pdf_plots:
    run: ../tools/files-to-folder.cwl
    in:
      input_files:
        source:
        - sc_rna_reduce/elbow_plot_pdf
        - sc_rna_reduce/qc_dim_corr_plot_pdf
        - sc_rna_reduce/umap_qc_mtrcs_plot_pdf
        - sc_rna_reduce/umap_plot_pdf
        - sc_rna_reduce/umap_spl_ph_plot_pdf
        - sc_rna_reduce/ccpca_plot_pdf
        - sc_rna_reduce/umap_spl_mito_plot_pdf
        - sc_rna_reduce/umap_spl_umi_plot_pdf
        - sc_rna_reduce/umap_spl_gene_plot_pdf
        - sc_rna_reduce/umap_spl_idnt_plot_pdf
        - sc_rna_reduce/ccpca_spl_idnt_plot_pdf
        - sc_rna_reduce/umap_spl_cnd_plot_pdf
        - sc_rna_reduce/umap_gr_cnd_spl_ph_plot_pdf
        - sc_rna_reduce/ccpca_spl_cnd_plot_pdf
        - sc_rna_reduce/umap_gr_cnd_spl_mito_plot_pdf
        - sc_rna_reduce/umap_gr_cnd_spl_umi_plot_pdf
        - sc_rna_reduce/umap_gr_cnd_spl_gene_plot_pdf
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

label: "Single-Cell RNA-Seq Dimensionality Reduction Analysis"
s:name: "Single-Cell RNA-Seq Dimensionality Reduction Analysis"
s:alternateName: "Integrates multiple single-cell RNA-Seq datasets, reduces dimensionality using PCA"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-rna-reduce.cwl
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
  Single-Cell RNA-Seq Dimensionality Reduction Analysis

  Removes noise and confounding sources of variation by reducing
  dimensionality of gene expression data from the outputs of
  “Single-Cell RNA-Seq Filtering Analysis” or “Single-Cell Multiome
  ATAC and RNA-Seq Filtering Analysis” pipelines. The results of
  this workflow are primarily used in “Single-Cell RNA-Seq Cluster
  Analysis” or “Single-Cell WNN Cluster Analysis” pipelines.