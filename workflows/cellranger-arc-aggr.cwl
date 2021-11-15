cwlVersion: v1.0
class: Workflow


requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
- class: MultipleInputFeatureRequirement


'sd:upstream':
  sc_rnaseq_sample:
  - "cellranger-arc-count.cwl"
  genome_indices:
    - "cellranger-mkref.cwl"

inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  gem_well_labels:
    type: string[]
    label: "scRNA-Seq Cell Ranger ARC Experiment"
    doc: "Array of GEM well identifiers to be used for labeling purposes only"
    'sd:upstreamSource': "sc_rnaseq_sample/alias"
    'sd:localLabel': true

  gex_molecule_info_h5:
    type: File[]
    label: "scRNA-Seq Cell Ranger ARC Experiment"
    doc: "Molecule-level information from individual runs of cellranger-arc count"
    'sd:upstreamSource': "sc_rnaseq_sample/gex_molecule_info_h5"

  atac_fragments_file:
    type: File[]
    secondaryFiles:
    - .tbi
    label: "scRNA-Seq Cell Ranger ARC Experiment"
    doc: "Count and barcode information from individual runs of cellranger-arc count"
    'sd:upstreamSource': "sc_rnaseq_sample/atac_fragments_file"

  barcode_metrics_report:
    type: File[]
    label: "scRNA-Seq Cell Ranger ARC Experiment"
    doc: "ATAC and GEX barcode metrics from individual runs of cellranger-arc count"
    'sd:upstreamSource': "sc_rnaseq_sample/barcode_metrics_report"

  indices_folder:
    type: Directory
    label: "Genome Type"
    doc: "Cell Ranger ARC generated genome indices folder"
    'sd:upstreamSource': "genome_indices/arc_indices_folder"

  normalization_mode:
    type:
    - "null"
    - type: enum
      symbols: ["none", "depth"]
    default: "none"
    label: "Library depth normalization mode"
    doc: "Library depth normalization mode"
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 4
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    'sd:layout':
      advanced: true

  memory_limit:
    type: int?
    default: 20
    label: "Maximum memory used (GB)"
    doc: "Maximum memory used (GB). The same will be applied to virtual memory"
    'sd:layout':
      advanced: true


outputs:

  web_summary_report:
    type: File
    outputSource: aggregate_counts/web_summary_report
    label: "Aggregated run summary metrics and charts in HTML format"
    doc: |
      Aggregated run summary metrics and charts in HTML format
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  metrics_summary_report:
    type: File
    outputSource: aggregate_counts/metrics_summary_report
    label: "Aggregated run summary metrics in CSV format"
    doc: |
      Aggregated run summary metrics in CSV format

  atac_fragments_file:
    type: File
    outputSource: aggregate_counts/atac_fragments_file
    label: "Aggregated count and barcode information"
    doc: |
      Count and barcode information for every ATAC fragment observed in the
      aggregated experiment in TSV format

  atac_peaks_bed_file:
    type: File
    outputSource: aggregate_counts/atac_peaks_bed_file
    label: "Locations of open-chromatin regions identified in aggregated experiment"
    doc: |
      Locations of open-chromatin regions identified in aggregated experiment
      (these regions are referred to as "peaks")

  atac_peak_annotation_file:
    type: File
    outputSource: aggregate_counts/atac_peak_annotation_file
    label: "Annotations of peaks based on genomic proximity alone for aggregated experiment"
    doc: |
      Annotations of peaks based on genomic proximity alone (for aggregated
      experiment). Note that these are not functional annotations and they
      do not make use of linkage with GEX data.

  secondary_analysis_report_folder:
    type: File
    outputSource: compress_secondary_analysis_report_folder/compressed_folder
    label: "Compressed folder with aggregated secondary analysis results"
    doc: |
      Compressed folder with secondary analysis results including dimensionality reduction,
      cell clustering, and differential expression of aggregated results

  filtered_feature_bc_matrix_folder:
    type: File
    outputSource: compress_filtered_feature_bc_matrix_folder/compressed_folder
    label: "Compressed folder with aggregated filtered feature-barcode matrices"
    doc: |
      Compressed folder with aggregated filtered feature-barcode matrices containing only cellular barcodes in MEX format

  filtered_feature_bc_matrix_h5:
    type: File
    outputSource: aggregate_counts/filtered_feature_bc_matrix_h5
    label: "Aggregated filtered feature-barcode matrices in HDF5 format"
    doc: |
      Aggregated filtered feature-barcode matrices containing only cellular barcodes in HDF5 format
  
  raw_feature_bc_matrices_folder:
    type: File
    outputSource: compress_raw_feature_bc_matrices_folder/compressed_folder
    label: "Compressed folder with aggregated unfiltered feature-barcode matrices"
    doc: |
      Compressed folder with aggregated unfiltered feature-barcode matrices containing all barcodes in MEX format

  raw_feature_bc_matrices_h5:
    type: File
    outputSource: aggregate_counts/raw_feature_bc_matrices_h5
    label: "Aggregated unfiltered feature-barcode matrices in HDF5 format"
    doc: |
      Aggregated unfiltered feature-barcode matrices containing all barcodes in HDF5 format

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
    label: "stdout log generated by cellranger-arc aggr"
    doc: |
      stdout log generated by cellranger-arc aggr

  aggregate_counts_stderr_log:
    type: File
    outputSource: aggregate_counts/stderr_log
    label: "stderr log generated by cellranger-arc aggr"
    doc: |
      stderr log generated by cellranger-arc aggr


steps:

  aggregate_counts:
    run: ../tools/cellranger-arc-aggr.cwl
    in:
      atac_fragments_file: atac_fragments_file
      barcode_metrics_report: barcode_metrics_report
      gex_molecule_info_h5: gex_molecule_info_h5
      gem_well_labels: gem_well_labels
      indices_folder: indices_folder
      normalization_mode: normalization_mode
      threads: threads
      memory_limit: memory_limit
      virt_memory_limit: memory_limit
    out:
    - web_summary_report
    - metrics_summary_report
    - atac_fragments_file
    - atac_peaks_bed_file
    - atac_peak_annotation_file
    - secondary_analysis_report_folder
    - filtered_feature_bc_matrix_folder
    - filtered_feature_bc_matrix_h5
    - raw_feature_bc_matrices_folder
    - raw_feature_bc_matrices_h5
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

  compress_raw_feature_bc_matrices_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: aggregate_counts/raw_feature_bc_matrices_folder
    out:
    - compressed_folder

  compress_secondary_analysis_report_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: aggregate_counts/secondary_analysis_report_folder
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Cell Ranger ARC Aggregate"
s:name: "Cell Ranger ARC Aggregate"
s:alternateName: "Aggregates data from multiple Cell Ranger ARC Count Gene Expression + ATAC experiments"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/cellranger-aggr.cwl
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
  Cell Ranger ARC Aggregate
  =========================
