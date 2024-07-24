cwlVersion: v1.0
class: Workflow


requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
- class: MultipleInputFeatureRequirement


"sd:upstream":
  sc_atacseq_sample:
  - "cellranger-atac-count.cwl"
  genome_indices:
  - "cellranger-mkref.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  gem_well_labels:
    type: string[]
    label: "Cell Ranger ATAC Count Experiment"
    doc: |
      Array of GEM well identifiers to be used for labeling purposes only.
      If not provided use rootnames of files from the barcode_metrics_report
      input
    "sd:upstreamSource": "sc_atacseq_sample/alias"
    "sd:localLabel": true

  fragments_file_from_count:
    type: File[]
    secondaryFiles:
    - .tbi
    label: "Cell Ranger ATAC Count Experiment"
    doc: |
      Array of files containing count and barcode information for
      every ATAC fragment observed in the "cellranger-atac count"
      experiment in TSV format.
    "sd:upstreamSource": "sc_atacseq_sample/atac_fragments_file"

  barcode_metrics_report_from_count:
    type: File[]
    label: "Cell Ranger ATAC Count Experiment"
    doc: |
      Array of files with per-barcode fragment counts & metrics
      produced by "cellranger-atac count" command in CSV format
    "sd:upstreamSource": "sc_atacseq_sample/barcode_metrics_report"

  indices_folder:
    type: Directory
    label: "Cell Ranger Reference Sample"
    doc: |
      Any "Cell Ranger Reference Sample" that
      builds a reference genome package of a
      selected species for quantifying gene
      expression and chromatin accessibility.
      This sample can be obtained from "Cell
      Ranger Reference (RNA, ATAC, RNA+ATAC)"
      pipeline.
    "sd:upstreamSource": "genome_indices/arc_indices_folder"
    "sd:localLabel": true

  annotation_gtf_file:
    type: File
    "sd:upstreamSource": "genome_indices/genome_indices/annotation_gtf"

  memory_limit:
    type: int?
    default: 20
    "sd:upstreamSource": "genome_indices/memory_limit"

  normalization_mode:
    type:
    - "null"
    - type: enum
      symbols: ["none", "depth"]
    default: "none"
    label: "Library depth normalization mode"
    doc: "Library depth normalization mode"
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

  web_summary_report:
    type: File
    outputSource: aggregate_counts/web_summary_report
    label: "Run summary metrics and charts in HTML format"
    doc: |
      Run summary metrics and charts in HTML format
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  metrics_summary_report_json:
    type: File
    outputSource: aggregate_counts/metrics_summary_report_json
    label: "Run summary metrics in JSON format"
    doc: |
      Run summary metrics in JSON format

  metrics_summary_report_csv:
    type: File
    outputSource: aggregate_counts/metrics_summary_report_csv
    label: "Run summary metrics in CSV format"
    doc: |
      Run summary metrics in CSV format

  barcode_metrics_report:
    type: File
    outputSource: aggregate_counts/barcode_metrics_report
    label: "Per-barcode fragment counts & metrics in CSV format"
    doc: |
      Per-barcode fragment counts & metrics in CSV format

  atac_fragments_file:
    type: File
    outputSource: aggregate_counts/fragments_file
    label: "Count and barcode information for every ATAC fragment"
    doc: |
      Count and barcode information for every ATAC fragment observed
      in the aggregated experiment in TSV format

  peaks_bed_file:
    type: File
    outputSource: aggregate_counts/peaks_bed_file
    label: "Locations of open-chromatin regions identified in the aggregated experiment"
    doc: |
      Locations of open-chromatin regions identified in the
      aggregated experiment (these regions are referred to
      as "peaks")

  peak_annotation_file:
    type: File
    outputSource: aggregate_counts/peak_annotation_file
    label: "Annotations of peaks based on genomic proximity alone"
    doc: |
      Annotations of peaks based on genomic proximity alone

  secondary_analysis_report_folder:
    type: File
    outputSource: compress_secondary_analysis_report_folder/compressed_folder
    label: "Folder with secondary analysis results"
    doc: |
      Folder with secondary analysis results

  filtered_feature_bc_matrix_folder:
    type: File
    outputSource: compress_filtered_feature_bc_matrix_folder/compressed_folder
    label: "Compressed folder with aggregated filtered peak-barcode matrices"
    doc: |
      Folder with aggregated filtered peak-barcode matrices
      containing only cellular barcodes in MEX format.

  filtered_feature_bc_matrix_h5:
    type: File
    outputSource: aggregate_counts/filtered_feature_bc_matrix_h5
    label: "Aggregated filtered peak-barcode matrices in HDF5 format"
    doc: |
      Aggregated filtered peak-barcode matrices containing
      only cellular barcodes in HDF5 format.

  # filtered_tf_bc_matrix_folder:
  #   type: File
  #   outputSource: compress_filtered_tf_bc_matrix_folder/compressed_folder
  #   label: "Compressed folder with aggregated filtered tf-barcode matrices"
  #   doc: |
  #     Folder with aggregated filtered tf-barcode matrices
  #     containing only cellular barcodes in MEX format.

  # filtered_tf_bc_matrix_h5:
  #   type: File
  #   outputSource: aggregate_counts/filtered_tf_bc_matrix_h5
  #   label: "Aggregated filtered tf-barcode matrices in HDF5 format"
  #   doc: |
  #     Aggregated filtered tf-barcode matrices containing
  #     only cellular barcodes in HDF5 format.

  loupe_browser_track:
    type: File
    outputSource: aggregate_counts/loupe_browser_track
    label: "Loupe Browser visualization and analysis file for aggregated results"
    doc: |
      Loupe Browser visualization and analysis file for aggregated results

  aggregation_metadata:
    type: File
    outputSource: aggregate_counts/aggregation_metadata
    label: "Aggregation metadata in CSV format"
    doc: |
      Aggregation metadata in CSV format

  aggregate_counts_stdout_log:
    type: File
    outputSource: aggregate_counts/stdout_log
    label: "stdout log generated by cellranger-atac aggr"
    doc: |
      stdout log generated by cellranger-atac aggr

  aggregate_counts_stderr_log:
    type: File
    outputSource: aggregate_counts/stderr_log
    label: "stderr log generated by cellranger-atac aggr"
    doc: |
      stderr log generated by cellranger-atac aggr

  html_data_folder:
    type: Directory
    outputSource: cellbrowser_build/html_data
    label: "Folder with not compressed CellBrowser formatted results"
    doc: |
      Folder with not compressed CellBrowser formatted results

  cellbrowser_report:
    type: File
    outputSource: cellbrowser_build/index_html_file
    label: "CellBrowser formatted Cellranger report"
    doc: |
      CellBrowser formatted Cellranger report
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"


steps:

  aggregate_counts:
    run: ../tools/cellranger-atac-aggr.cwl
    in:
      fragments_file_from_count: fragments_file_from_count
      barcode_metrics_report: barcode_metrics_report_from_count
      gem_well_labels: gem_well_labels
      indices_folder: indices_folder
      normalization_mode: normalization_mode
      threads:
        source: threads
        valueFrom: $(parseInt(self))
      memory_limit: memory_limit
      virt_memory_limit: memory_limit
    out:
    - web_summary_report
    - metrics_summary_report_json
    - metrics_summary_report_csv
    - barcode_metrics_report
    - fragments_file
    - peaks_bed_file
    - peak_annotation_file
    - secondary_analysis_report_folder
    - filtered_feature_bc_matrix_folder
    - filtered_feature_bc_matrix_h5
    # - filtered_tf_bc_matrix_folder
    # - filtered_tf_bc_matrix_h5
    - aggregation_metadata
    - loupe_browser_track
    - stdout_log
    - stderr_log

  compress_filtered_feature_bc_matrix_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: aggregate_counts/filtered_feature_bc_matrix_folder
    out:
    - compressed_folder

  # compress_filtered_tf_bc_matrix_folder:
  #   run: ../tools/tar-compress.cwl
  #   in:
  #     folder_to_compress: aggregate_counts/filtered_tf_bc_matrix_folder
  #   out:
  #   - compressed_folder

  compress_secondary_analysis_report_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: aggregate_counts/secondary_analysis_report_folder
    out:
    - compressed_folder

  cellbrowser_build:
    run: ../tools/cellbrowser-build-cellranger-atac.cwl
    in:
      secondary_analysis_report_folder: aggregate_counts/secondary_analysis_report_folder
      filtered_feature_bc_matrix_folder: aggregate_counts/filtered_feature_bc_matrix_folder
      aggregation_metadata: aggregate_counts/aggregation_metadata
      annotation_gtf_file: annotation_gtf_file
    out:
    - html_data
    - index_html_file

$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Cell Ranger Aggregate (ATAC)"
s:name: "Cell Ranger Aggregate (ATAC)"
s:alternateName: "Combines outputs from multiple runs of Cell Ranger Count (ATAC) pipeline"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/cellranger-atac-aggr.cwl
s:codeRepository: https://github.com/datirium/workflows
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
  Cell Ranger Aggregate (ATAC)

  Combines outputs from multiple runs of
  “Cell Ranger Count (ATAC)” pipeline.