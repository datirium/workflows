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
    - var split_by_comma = function(line) {
          function get_unique(value, index, self) {
            return self.indexOf(value) === index && value != "";
          }
          var splitted_line = line?line.split(/,+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };


"sd:upstream":
  sc_tools_sample:
  - "sc-atac-reduce.cwl"
  - "sc-atac-cluster.cwl"
  - "sc-rna-reduce.cwl"
  - "sc-rna-cluster.cwl"
  - "sc-wnn-cluster.cwl"
  - "sc-ctype-assign.cwl"
  - "sc-rna-azimuth.cwl"
  sc_atac_sample:
  - "cellranger-arc-count.cwl"
  - "cellranger-arc-aggr.cwl"
  - "cellranger-atac-count.cwl"
  - "cellranger-atac-aggr.cwl"
  genome_indices:
  - "genome-indices.cwl"


inputs:

  alias:
    type: string
    label: "Analysis name"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Single-cell Analysis with ATAC-Seq Datasets"
    doc: |
      Analysis that includes single-cell
      multiome RNA and ATAC-Seq or just
      ATAC-Seq datasets run through either
      "Single-Cell ATAC-Seq Dimensionality
      Reduction Analysis", "Single-Cell
      RNA-Seq Dimensionality Reduction
      Analysis", "Single-Cell WNN Cluster
      Analysis" or "Single-Cell RNA-Seq
      Reference Mapping" pipeline at any
      of the processing stages.
    "sd:upstreamSource": "sc_tools_sample/seurat_data_rds"
    "sd:localLabel": true

  query_reduction:
    type:
    - "null"
    - type: enum
      symbols:
      - "RNA"
      - "ATAC"
      - "WNN"
      - "REF"                                # from the sc-rna-azimuth.cwl pipeline
    default: "RNA"
    label: "Dimensionality reduction"
    doc: |
      Dimensionality reduction to be used
      for generating UMAP plots.

  atac_fragments_file:
    type: File
    secondaryFiles:
    - .tbi
    label: "Cell Ranger ATAC or RNA+ATAC Sample"
    doc: |
      Any "Cell Ranger ATAC or RNA+ATAC Sample"
      for loading chromatin accessibility data
      from. This sample can be obtained from
      one of the following pipelines: "Cell
      Ranger Count (RNA+ATAC)", "Cell Ranger
      Aggregate (RNA+ATAC)", "Cell Ranger Count
      (ATAC)", or "Cell Ranger Aggregate (ATAC)".
    "sd:upstreamSource": "sc_atac_sample/atac_fragments_file"
    "sd:localLabel": true

  genome_type:
    type: string
    label: "Genome type"
    doc: |
      Reference genome
    "sd:upstreamSource": "genome_indices/genome"
    "sd:localLabel": true

  annotation_file:
    type: File
    "sd:upstreamSource": "genome_indices/annotation"

  chrom_length_file:                                            # not used, but removing it prevents SciDAP from rerunning old samples. Keeping it for now.
    type: File?
    "sd:upstreamSource": "genome_indices/chrom_length"

  blacklist_regions_file:
    type:
      type: enum
      symbols:
      - "hg19"
      - "hg38"
      - "mm10"
    "sd:upstreamSource": "genome_indices/genome"

  groupby:
    type: string?
    default: null
    label: "Subsetting category (optional)"
    doc: |
      Single cell metadata column to group
      cells for optional subsetting before
      running differential accessibility analysis.
      To group cells by dataset, use "dataset".
      Custom groups can be defined based on
      any single cell metadata added through
      the "Datasets metadata (optional)" or
      "Selected cell barcodes (optional)"
      inputs. Default: do not subset cells

  subset:
    type: string?
    default: null
    label: "Subsetting values (optional)"
    doc: |
      Comma separated list of values from
      the single cell metadata column
      selected in "Subsetting category
      (optional)" input. Ignored if grouping
      category is not provided. Default: do
      not subset cells

  splitby:
    type: string
    label: "Comparison category"
    doc: |
      Single cell metadata column to split
      cells into two comparison groups before
      running differential accessibility analysis.
      To split cells by dataset, use "dataset".
      Custom groups can be defined based on
      any single cell metadata added through
      the "Datasets metadata (optional)" or
      "Selected cell barcodes (optional)"
      inputs. The direction of comparison is
      always "Second comparison group" vs
      "First comparison group".

  first_cond:
    type: string
    label: "First comparison group"
    doc: |
      Value from the single cell metadata
      column selected in "Comparison category"
      input to define the first group of cells
      for differential accessibility analysis.

  second_cond:
    type: string
    label: "Second comparison group"
    doc: |
      Value from the single cell metadata
      column selected in "Comparison category"
      input to define the second group of cells
      for differential accessibility analysis.

  analysis_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "negative-binomial"                              # (negbinom) Negative Binomial Generalized Linear Model (use FindMarkers with peaks from Seurat object)
      - "poisson"                                        # (poisson) Poisson Generalized Linear Model (use FindMarkers with peaks from Seurat object)
      - "logistic-regression"                            # (LR) Logistic Regression (use FindMarkers with peaks from Seurat object)
      - "mast"                                           # (MAST) MAST package (use FindMarkers with peaks from Seurat object)
      - "manorm2-peaks-by-dataset"                       # call peaks for each dataset with MACS2, then run MAnorm2 with datasets
      - "manorm2-peaks-by-comparison-category"           # call peaks for each comparison group with MACS2, then run MAnorm2 with datasets
      - "diffbind-deseq-peaks-by-dataset"                # call peaks for each dataset with MACS2, then run DiffBind (DESeq2) with datasets
      - "diffbind-deseq-peaks-by-comparison-category"    # call peaks for each comparison group with MACS2, then run DiffBind (DESeq2) with datasets
      - "diffbind-edger-peaks-by-dataset"                # call peaks for each dataset with MACS2, then run DiffBind (EdgeR) with datasets
      - "diffbind-edger-peaks-by-comparison-category"    # call peaks for each comparison group with MACS2, then run DiffBind (EdgeR) with datasets
    default: "logistic-regression"
    label: "Statistical test"
    doc: |
      Statistical test to use in the
      differential accessibility analysis.
      Chromatin accessibility data will first
      be filtered to include fragments from the
      cells retained in the loaded single-cell
      datasets after optional filtering by the
      "Subsetting category/values" and/or
      "Selected cell barcodes". The resulting
      fragments will then be splitted by the
      "Comparison category" and, if
      "*-peaks-by-dataset" or
      "*-peaks-by-comparison-category"
      is selected, they will be aggregated to
      the pseudo bulk form for peak calling
      with MACS2, either by dataset or by
      comparison category, respectively.
      Othwerwise, the analysis will be run on
      the cells level using peaks already
      available in the loaded single-cell
      datasets.
      Default: logistic-regression

  maximum_padj:
    type: float?
    default: 0.05
    label: "Maximum adjusted p-value"
    doc: |
      Maximum adjusted p-value threshold
      for selecting differentially accessible
      regions to be visualized on the heatmap.
      Default: 0.05

  minimum_logfc:
    type: float?
    default: 1
    label: "Minimum log2 fold change absolute value"
    doc: |
      Minimum log2 fold change absolute value
      for selecting differentially accessible
      regions to be visualized on the heatmap.
      Default: 1.0

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
      "datasets_metadata.tsv" output from
      the "Files" tab of either "Single-Cell
      Multiome ATAC-Seq and RNA-Seq Filtering
      Analysis" or "Single-Cell ATAC-Seq
      Filtering Analysis" and add extra
      columns as needed. The pipeline selection
      depends on whether the "Single-cell
      Analysis with ATAC-Seq Datasets" includes
      both RNA and ATAC-Seq or only ATAC-Seq
      datasets.

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
      from "Single-cell Analysis with ATAC-Seq
      Datasets" and can be utilized in the
      current or future steps of analysis.

  minimum_qvalue:
    type: float?
    default: 0.05
    label: "Minimum MACS2 peak calling FDR (ignored if statistical test is not using MAnorm2 or DiffBind)"
    doc: |
      Minimum FDR (q-value) cutoff for MACS2
      peak detection. Ignored if "Statistical
      test" input is not set to "manorm2-*"
      or "diffbind-*" methods.
      Default: 0.05
    "sd:layout":
      advanced: true

  maximum_peaks:
    type: int?
    default: 0
    label: "Maximum number of peaks to keep, set to 0 to keep all (ignored if statistical test is not using MAnorm2)"
    doc: |
      The maximum number of the most significant
      peaks to select for each of the comparison
      groups when constructing reference genomic
      bins. The top significant peaks are selected
      based on the score column which is calculated
      by MACS2 as int(-10*log10qvalue). Ignored
      if "Statistical test" input is not set to
      "manorm2-peaks-by-dataset" or
      "manorm2-peaks-by-comparison-category".
      Default: 0 - keep all peaks.
    "sd:layout":
      advanced: true

  minimum_peak_gap:
    type: int?
    default: 150
    label: "Minimum distance for merging peaks (ignored if statistical test is not using MAnorm2)"
    doc: |
      Peaks remained after optional filtering
      by the "Maximum number of peaks to keep"
      will be merged if the distance between
      them is smaller than the provided value.
      Merging is first applied within comparison
      group and then to all peaks together before
      splitting them into the reference genomic
      bins of the selected size. Ignored if
      "Statistical test" input is not set to
      "manorm2-peaks-by-dataset" or
      "manorm2-peaks-by-comparison-category".
      Default: 150
    "sd:layout":
      advanced: true

  bin_size:
    type: int?
    default: 1000
    label: "Bin size (ignored if statistical test is not using MAnorm2)"
    doc: |
      The size of non-overlapping reference
      genomic bins used by MAnorm2 when
      generating a table of reads per bins
      counts. Ignored if "Statistical test"
      input is not set to
      "manorm2-peaks-by-dataset" or
      "manorm2-peaks-by-comparison-category".
      Default: 1000
    "sd:layout":
      advanced: true

  minimum_overlap:
    type: float?
    default: 0.5
    label: "Minimum overlap fraction between the datasets (ignored if statistical test is not using MAnorm2 or DiffBind)"
    doc: |
      Keep only those reference genomic
      bins that are present in at least
      this fraction of datasets within
      each of the comparison groups.
      Ignored if "Statistical test"
      input is not set to "manorm2-*"
      or "diffbind-*" methods.
      Default: 0.5
    "sd:layout":
      advanced: true

  promoter_dist:
    type: int?
    default: 1000
    label: "Promoter distance, bp"
    doc: |
      Max distance from the gene TSS (in both
      direction) overlapping which the
      differentially accessible site will be
      assigned to the promoter region.
      Default: 1000 bp
    "sd:layout":
      advanced: true

  upstream_dist:
    type: int?
    default: 20000
    label: "Upstream distance, bp"
    doc: |
      Max distance from the promoter (only
      in upstream direction) overlapping
      which the differentially accessible site
      will be assigned to the upstream region.
      Default: 20,000 bp
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
    label: "Number of cores/cpus to use"
    doc: |
      Parallelization parameter to define the
      number of cores/CPUs that can be utilized
      simultaneously.
      Default: 4
    "sd:layout":
      advanced: true


outputs:

  all_db_sites_tsv_with_labels:
    type: File
    outputSource: add_label_column/output_file
    label: "Differentially accessible regions with labels (not filtered)"
    doc: |
      Not filtered by "Maximum adjusted p-value"
      and "Minimum log2 fold change absolute value"
      differentially accessible regions with labels.
      TSV format.
    "sd:visualPlugins":
    - queryRedirect:
        tab: "Overview"
        label: "Volcano Plot"
        url: "https://scidap.com/vp/volcano"
        query_eval_string: "`data_file=${this.getSampleValue('outputs', 'all_db_sites_tsv_with_labels')}&data_col=label&x_col=log2FoldChange&y_col=padj`"

  tag_dnst_htmp_html:
    type: File?
    outputSource: extend_htmp/tag_dnst_htmp_html
    label: "Tag Density Heatmap"
    doc: |
      Tag density heatmap around the centers
      of differentially accessible regions.
      Filtered by "Maximum adjusted p-value"
      and "Minimum log2 fold change absolute
      value"; optionally subsetted to include
      only cells with "Subsetting values
      (optional)" from the "Subsetting category
      (optional)".
      HTML format.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  sc_report_html_file:
    type: File?
    outputSource: sc_atac_dbinding/sc_report_html_file
    label: "Analysis log"
    doc: |
      Tehcnical report.
      HTML format.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  experiment_info:
    type: File
    label: "Genome coverage tracks order"
    doc: |
      Markdown file to explain
      the order of genome coverage
      tracks in the IGV browser.
    outputSource: create_metadata/output_file
    "sd:visualPlugins":
    - markdownView:
        tab: "Overview"

  cell_cnts_plot_png:
    type: File?
    outputSource: sc_atac_dbinding/cell_cnts_plot_png
    label: "Number of cells per dataset or comparison group"
    doc: |
      Number of cells per dataset or comparison
      group. Colored by "Comparison category";
      optionally subsetted to include only cells
      with "Subsetting values (optional)" from
      the "Subsetting category (optional)".
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Number of cells per dataset or comparison group"

  umap_spl_tst_plot_png:
    type: File?
    outputSource: sc_atac_dbinding/umap_spl_tst_plot_png
    label: "UMAP colored by selected for analysis cells"
    doc: |
      UMAP colored by selected for analysis cells.
      Split by the single cell metadata column
      defined in the "Comparison category";
      optionally subsetted to include only cells
      with "Subsetting values (optional)" from
      the "Subsetting category (optional)".
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "UMAP colored by selected for analysis cells"

  vlcn_plot_png:
    type: File?
    outputSource: sc_atac_dbinding/vlcn_plot_png
    label: "Volcano plot of differentially accessible regions"
    doc: |
      Volcano plot of differentially
      accessible regions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Volcano plot"
        Caption: "Volcano plot of differentially accessible regions"

  tag_dnst_htmp_plot_png:
    type: File?
    outputSource: sc_atac_dbinding/tag_dnst_htmp_plot_png
    label: "Tag density heatmap around the centers of filtered diff. accessible regions"
    doc: |
      Tag density heatmap around the centers
      of differentially accessible regions.
      Filtered by "Maximum adjusted p-value"
      and "Minimum log2 fold change absolute
      value"; optionally subsetted to include
      only cells with "Subsetting values
      (optional)" from the "Subsetting category
      (optional)".
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Heatmap"
        Caption: "Tag density heatmap around the centers of filtered diff. accessible regions"

  tag_dnst_htmp_gct:
    type: File?
    outputSource: extend_htmp/tag_dnst_htmp_gct
    label: "Tag density heatmap around the centers of filtered diff. accessible regions"
    doc: |
      Tag density heatmap around the centers
      of differentially accessible regions.
      Filtered by "Maximum adjusted p-value"
      and "Minimum log2 fold change absolute
      value"; optionally subsetted to include
      only cells with "Subsetting values
      (optional)" from the "Subsetting category
      (optional)".
      GCT format.

  tag_dnst_htmp_tsv:
    type: File?
    outputSource: extend_htmp/tag_dnst_htmp_tsv
    label: "Tag density heatmap around the centers of filtered diff. accessible regions"
    doc: |
      Tag density heatmap around the centers
      of differentially accessible regions.
      Filtered by "Maximum adjusted p-value"
      and "Minimum log2 fold change absolute
      value"; optionally subsetted to include
      only cells with "Subsetting values
      (optional)" from the "Subsetting category
      (optional)".
      TSV format.

  coverage_first_bigwig_file:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_dbinding/coverage_first_bigwig_file
    label: "Genome coverage (first)"
    doc: |
      Normalized genome coverage calculated
      from either ATAC fragments or extended
      to 40bp lenght Tn5 cut sites. The data
      were first optionally subsetted to include
      only cells with "Subsetting values (optional)"
      from the "Subsetting category (optional)",
      then split either by dataset or by the
      "Comparison category" and aggregated to
      the pseudobulk form. First comparison group.
      BigWig format.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "wig"
        name: "Genome coverage (first)"
        height: 120

  coverage_second_bigwig_file:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_dbinding/coverage_second_bigwig_file
    label: "Genome coverage (second)"
    doc: |
      Normalized genome coverage calculated
      from either ATAC fragments or extended
      to 40bp lenght Tn5 cut sites. The data
      were first optionally subsetted to include
      only cells with "Subsetting values (optional)"
      from the "Subsetting category (optional)",
      then split either by dataset or by the
      "Comparison category" and aggregated to
      the pseudobulk form. Second comparison group.
      BigWig format.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "wig"
        name: "Genome coverage (second)"
        height: 120

  peaks_first_bed_file:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_dbinding/peaks_first_bed_file
    label: "Called peaks (first)"
    doc: |
      Peaks called by MACS2 from the Tn5 cut
      sites split either by dataset or by
      "Comparison category". Optionally
      subsetted to include only cells with
      "Subsetting values (optional)" from
      the "Subsetting category (optional)".
      First comparison group.
      NarrowPeak format.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        name: "Called peaks (first)"
        displayMode: "COLLAPSE"
        height: 40

  peaks_second_bed_file:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_dbinding/peaks_second_bed_file
    label: "Called peaks (second)"
    doc: |
      Peaks called by MACS2 from the Tn5 cut
      sites split either by dataset or by
      "Comparison category". Optionally
      subsetted to include only cells with
      "Subsetting values (optional)" from
      the "Subsetting category (optional)".
      Second comparison group.
      NarrowPeak format.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        name: "Called peaks (second)"
        displayMode: "COLLAPSE"
        height: 40

  peaks_first_xls_file:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_dbinding/peaks_first_xls_file
    label: "Called peaks (first, xls format)"
    doc: |
      Peaks called by MACS2 from the Tn5 cut
      sites split either by dataset or by
      "Comparison category". Optionally
      subsetted to include only cells with
      "Subsetting values (optional)" from
      the "Subsetting category (optional)".
      First comparison group.
      XLS format.

  peaks_second_xls_file:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_dbinding/peaks_second_xls_file
    label: "Called peaks (second, xls format)"
    doc: |
      Peaks called by MACS2 from the Tn5 cut
      sites split either by dataset or by
      "Comparison category". Optionally
      subsetted to include only cells with
      "Subsetting values (optional)" from
      the "Subsetting category (optional)".
      Second comparison group.
      XLS format.

  summits_first_bed_file:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_dbinding/summits_first_bed_file
    label: "Peaks summits (first)"
    doc: |
      Summits of the peaks called by MACS2
      from the Tn5 cut sites split either
      by dataset or by "Comparison category".
      Optionally subsetted to include only
      cells with "Subsetting values (optional)"
      from the "Subsetting category (optional)".
      First comparison group.
      BED format.

  summits_second_bed_file:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_atac_dbinding/summits_second_bed_file
    label: "Peaks summits (second)"
    doc: |
      Summits of the peaks called by MACS2
      from the Tn5 cut sites split either
      by dataset or by "Comparison category".
      Optionally subsetted to include only
      cells with "Subsetting values (optional)"
      from the "Subsetting category (optional)".
      Second comparison group.
      BED format.

  dflt_peaks_bigbed_file:
    type: File?
    outputSource: sc_atac_dbinding/dflt_peaks_bigbed_file
    label: "Seurat or Cell Ranger peaks"
    doc: |
      Peaks extracted from the loaded
      "Single-cell Analysis with ATAC-Seq
      Datasets".
      BigBed format.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        format: "bigbed"
        name: "Seurat or Cell Ranger peaks"
        displayMode: "COLLAPSE"
        height: 40

  all_db_sites_tsv:
    type: File
    outputSource: restore_columns/output_file
    label: "Differentially accessible regions (not filtered)"
    doc: |
      Not filtered by "Maximum adjusted p-value"
      and "Minimum log2 fold change absolute value"
      differentially accessible regions with the
      nearest genes assigned.
      TSV format.
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Diff. accessible regions"
        Title: "Differentially accessible regions (not filtered)"

  all_db_sites_bed:
    type: File
    outputSource: sc_atac_dbinding/all_db_sites_bed
    label: "Differentially accessible regions (not filtered)"
    doc: |
      Not filtered by "Maximum adjusted p-value"
      and "Minimum log2 fold change absolute value"
      differentially accessible regions.
      BED format.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        name: "Differentially accessible regions (not filtered)"
        displayMode: "COLLAPSE"
        height: 40

  fltr_db_sites_bed:
    type: File?
    outputSource: sc_atac_dbinding/fltr_db_sites_bed
    label: "Differentially accessible regions (filtered)"
    doc: |
      Filtered by "Maximum adjusted p-value"
      and "Minimum log2 fold change absolute
      value" differentially accessible regions.
      BED format.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        name: "Differentially accessible regions (filtered)"
        displayMode: "COLLAPSE"
        height: 40

  pdf_plots:
    type: File
    outputSource: compress_pdf_plots/compressed_folder
    label: "Compressed folder with all PDF plots"
    doc: |
      Compressed folder with all PDF plots.

  sc_atac_dbinding_human_log:
    type: File?
    outputSource: sc_atac_dbinding/human_log
    label: "Human readable error log"
    doc: |
      Human readable error log
      from the sc_atac_dbinding step.

  sc_atac_dbinding_stdout_log:
    type: File
    outputSource: sc_atac_dbinding/stdout_log
    label: "Output log"
    doc: |
      Stdout log from the
      sc_atac_dbinding step.

  sc_atac_dbinding_stderr_log:
    type: File
    outputSource: sc_atac_dbinding/stderr_log
    label: "Error log"
    doc: |
      Stderr log from the
      sc_atac_dbinding step.


steps:

  sc_atac_dbinding:
    run: ../tools/sc-atac-dbinding.cwl
    in:
      query_data_rds: query_data_rds
      reduction:
        source: query_reduction
        valueFrom: |
          ${
            if (self == "RNA") {
              return "rnaumap";
            } else if (self == "ATAC") {
              return "atacumap";
            } else if (self == "WNN") {
              return "wnnumap";
            } else {
              return "refumap";
            }
          }
      atac_fragments_file: atac_fragments_file
      datasets_metadata: datasets_metadata
      barcodes_data: barcodes_data
      groupby:
        source: groupby
        valueFrom: |
          ${
            if (self == "dataset") {
              return "new.ident";
            } else if (self == "") {
              return null;
            } else {
              return self;
            }
          }
      subset:
        source: subset
        valueFrom: $(split_by_comma(self))
      splitby:
        source: splitby
        valueFrom: |
          ${
            if (self == "dataset") {
              return "new.ident";
            } else {
              return self;
            }
          }
      first_cond: first_cond
      second_cond: second_cond
      analysis_method:
        source: analysis_method
        valueFrom: |
          ${
            if (self == "manorm2-peaks-by-dataset") {
              return "manorm2-full";
            } else if (self == "manorm2-peaks-by-comparison-category") {
              return "manorm2-half";
            } else if (self == "diffbind-deseq-peaks-by-dataset") {
              return "diffbind-deseq-full";
            } else if (self == "diffbind-deseq-peaks-by-comparison-category") {
              return "diffbind-deseq-half";
            } else if (self == "diffbind-edger-peaks-by-dataset") {
              return "diffbind-edger-full";
            } else if (self == "diffbind-edger-peaks-by-comparison-category") {
              return "diffbind-edger-half";
            } else {
              return self;
            }
          }
      genome_type:
        source: genome_type
        valueFrom: $(self=="mm10"?"mm":"hs")
      minimum_qvalue: minimum_qvalue
      minimum_peak_gap: minimum_peak_gap
      bin_size: bin_size
      minimum_overlap: minimum_overlap
      maximum_peaks:
        source: maximum_peaks
        valueFrom: $(self==0?null:self)                 # to return null for 0
      blacklist_regions_file: blacklist_regions_file
      maximum_padj: maximum_padj
      minimum_logfc: minimum_logfc
      verbose:
        default: true
      export_pdf_plots:
        default: true
      parallel_memory_limit:
        default: 32
      vector_memory_limit:
        default: 128
      export_html_report: export_html_report
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - umap_spl_tst_plot_png
    - cell_cnts_plot_png
    - vlcn_plot_png
    - tag_dnst_htmp_plot_png
    - tag_dnst_htmp_gct
    - all_db_sites_tsv
    - all_db_sites_bed
    - fltr_db_sites_bed
    - coverage_first_bigwig_file
    - coverage_second_bigwig_file
    - peaks_first_bed_file
    - peaks_second_bed_file
    - peaks_first_xls_file
    - peaks_second_xls_file
    - summits_first_bed_file
    - summits_second_bed_file
    - dflt_peaks_bigbed_file
    - all_plots_pdf
    - sc_report_html_file
    - human_log
    - stdout_log
    - stderr_log

  folder_pdf_plots:
    run: ../tools/files-to-folder.cwl
    in:
      input_files:
        source:
        - sc_atac_dbinding/all_plots_pdf
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

  filter_columns:
    run: ../tools/custom-bash.cwl
    in:
      input_file: sc_atac_dbinding/all_db_sites_tsv
      script:
        default: >
          cat $0 | grep -v "start" | awk
          'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"}
          {print $1"\t"$2"\t"$3"\t"$3-$2+1"\t0\t"NR"\t0\t0\t0\t0"}' > `basename $0`
    out:
    - output_file

  assign_genes:
    run: ../tools/iaintersect.cwl
    in:
      input_filename: filter_columns/output_file
      annotation_filename: annotation_file
      promoter_bp: promoter_dist
      upstream_bp: upstream_dist
    out:
    - result_file

  restore_columns:
    run: ../tools/custom-bash.cwl
    in:
      input_file:
      - assign_genes/result_file
      - sc_atac_dbinding/all_db_sites_tsv
      script:
        default: |
          cat $0 | grep -v "start" | sort -k 11n | cut -f 1-5,15 > tmp_iaintersect_result.tsv
          cat $1 | grep -v "start" > tmp_sc_atac_dbinding_result.tsv
          HEADER=`head -n 1 $1`;
          echo -e "refseq_id\tgene_id\ttxStart\ttxEnd\tstrand\tregion\t${HEADER}" > `basename $1`;
          cat tmp_iaintersect_result.tsv | paste - tmp_sc_atac_dbinding_result.tsv >> `basename $1`
          rm tmp_iaintersect_result.tsv tmp_sc_atac_dbinding_result.tsv
    out:
    - output_file

  add_label_column:
    run: ../tools/custom-bash.cwl
    in:
      input_file: restore_columns/output_file
      script:
        default: |
          HEADER=`head -n 1 $0`;
          echo -e "label\t${HEADER}" > sc_all_db_sites_labeled.tsv;
          cat "$0" | grep -v "start" | awk -F "\t" '{print $7":"$8"-"$9" "$2" "$6"\t"$0}' >> sc_all_db_sites_labeled.tsv
    out:
    - output_file

  extend_htmp:
    run: ../tools/sc-utils-extend-htmp.cwl
    in:
      gct_file: sc_atac_dbinding/tag_dnst_htmp_gct
      tsv_file: add_label_column/output_file
    out:
    - tag_dnst_htmp_gct
    - tag_dnst_htmp_html
    - tag_dnst_htmp_tsv

  create_metadata:
    run: ../tools/custom-bash.cwl
    in:
      input_file:
        source:
        - sc_atac_dbinding/coverage_first_bigwig_file
        - sc_atac_dbinding/coverage_second_bigwig_file
        valueFrom: $(self.flat())
      param:
        source:
        - sc_atac_dbinding/coverage_first_bigwig_file
        - sc_atac_dbinding/coverage_second_bigwig_file
        valueFrom: $(self[0].length + "," + self[1].length)
      script:
        default: |
          #!/bin/bash
          set -- "$0" "$@"
          IFS=',' read -r FIRST_N SECOND_N <<< "${!#}"
          FILES=("${@:1:$(($#-1))}")
          echo "| File | Position |" > experiment_info.md
          echo "| :-- | --: |" >> experiment_info.md
          echo "| _First comparison group_ | |" >> experiment_info.md
          for ((i = 0; i < FIRST_N; i++)); do
            echo "| $(basename "${FILES[i]}") | $((i+1)) |" >> experiment_info.md
          done
          echo "| _Second comparison group_ | |" >> experiment_info.md
          for ((i = 0; i < SECOND_N; i++)); do
            j=$((FIRST_N + i))
            echo "| $(basename "${FILES[j]}") | $((i+1)) |" >> experiment_info.md
          done
    out:
    - output_file


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-Cell ATAC-Seq Differential Accessibility Analysis"
s:name: "Single-Cell ATAC-Seq Differential Accessibility Analysis"
s:alternateName: "Identifies differentially accessible regions between two groups of cells"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-atac-dbinding.cwl
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
  Single-Cell ATAC-Seq Differential Accessibility Analysis

  Identifies differentially accessible regions between any
  two groups of cells, optionally aggregating chromatin
  accessibility data from single-cell to pseudo bulk form.