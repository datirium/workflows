cwlVersion: v1.0
class: Workflow


requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
- class: MultipleInputFeatureRequirement


"sd:upstream":
  sc_experiment:
  - "single-cell-preprocess-cellranger.cwl"
  - "cellranger-multi.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  molecule_info_h5:
    type:
    - "null"
    -  File[]
    label: "Cell Ranger RNA or RNA+VDJ Sample"
    doc: |
      Any "Cell Ranger RNA or RNA+VDJ Sample"
      that produces gene expression and,
      optionally, V(D)J contigs data, from a
      single 10x Genomics library
    "sd:upstreamSource": "sc_experiment/molecule_info_h5"
    "sd:localLabel": true

  filtered_data_folder:
    type:
    - "null"
    - Directory[]
    "sd:upstreamSource": "sc_experiment/filtered_data_folder"

  gem_well_labels:
    type: string[]
    "sd:upstreamSource": "sc_experiment/alias"

  normalization_mode:
    type:
    - "null"
    - type: enum
      symbols:
      - "none"
      - "mapped"
    default: "none"
    label: "Library depth normalization mode"
    doc: "Library depth normalization mode"
    "sd:layout":
      advanced: true

  clonotype_grouping:
    type:
    - "null"
    - type: enum
      name: "clonotype_grouping"
      symbols:
      - "same_donor_different_origins"
      - "same_donor_and_origin"
      - "different_donors"
    default: "different_donors"
    label: "Clonotype grouping. Ignored if upstream analysis doesn't include V(D)J data"
    doc: |
      When cellranger aggr is called with cellranger multi outputs, there are three
      ways it can process the datasets depending on the combination of donor and
      origin values
    "sd:layout":
      advanced: true

  threads:
    type: int?
    default: 4
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    "sd:layout":
      advanced: true

  memory_limit:
    type: int?
    default: 30
    label: "Maximum memory used (GB)"
    doc: "Maximum memory used (GB). The same will be applied to virtual memory"
    "sd:layout":
      advanced: true


outputs:

  web_summary_report:
    type: File
    outputSource: aggregate_counts/web_summary_report
    label: "Aggregated run summary metrics and charts in HTML format"
    doc: "Aggregated run summary metrics and charts in HTML format"
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  metrics_summary_report_json:
    type: File
    outputSource: aggregate_counts/metrics_summary_report_json
    label: "Aggregated GEX run summary metrics in JSON format"
    doc: "Aggregated GEX run summary metrics in JSON format"

  secondary_analysis_report_folder:
    type: File
    outputSource: compress_secondary_analysis_report_folder/compressed_folder
    label: "Compressed folder with aggregated secondary analysis results"
    doc: |
      Compressed folder with secondary analysis of GEX data including dimensionality
      reduction, cell clustering, and differential expression

  filtered_feature_bc_matrix_folder:
    type: File
    outputSource: compress_filtered_feature_bc_matrix_folder/compressed_folder
    label: "Compressed folder with aggregated filtered feature-barcode matrices"
    doc: |
      Compressed folder with aggregated filtered feature-barcode matrices
      containing only cellular barcodes in MEX format

  filtered_feature_bc_matrix_h5:
    type: File
    outputSource: aggregate_counts/filtered_feature_bc_matrix_h5
    label: "Aggregated filtered feature-barcode matrices in HDF5 format"
    doc: |
      Filtered feature-barcode matrices containing only cellular
      barcodes in HDF5 format

  aggregation_metadata:
    type: File
    outputSource: aggregate_counts/aggregation_metadata
    label: "Aggregation metadata in CSV format"
    doc: "Aggregation metadata in CSV format"

  grouping_data:
    type: File
    outputSource: aggregate_counts/grouping_data
    label: "Example of datasets grouping"
    doc: "Example of TSV file to define datasets grouping"

  loupe_browser_track:
    type: File
    outputSource: aggregate_counts/loupe_browser_track
    label: "Loupe Browser visualization and analysis file"
    doc: "Loupe Browser visualization and analysis file"

  clonotypes_csv:
    type: File?
    outputSource: aggregate_counts/clonotypes_csv
    label: "CSV file with high-level descriptions of each clonotype"
    doc: "CSV file with high-level descriptions of each clonotype"
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "V(D)J clonotypes"
        Title: "V(D)J clonotypes"

  consensus_sequences_fasta:
    type: File?
    outputSource: aggregate_counts/consensus_sequences_fasta
    label: "The consensus sequence of each assembled contig"
    doc: "The consensus sequence of each assembled contig"

  consensus_annotations_csv:
    type: File?
    outputSource: aggregate_counts/consensus_annotations_csv
    label: "CSV file with high-level and detailed annotations of each clonotype consensus sequence"
    doc: "CSV file with high-level and detailed annotations of each clonotype consensus sequence"

  filtered_contig_annotations_csv:
    type: File?
    outputSource: aggregate_counts/filtered_contig_annotations_csv
    label: "CSV file with high-level annotations of each high-confidence contig from cell-associated barcodes"
    doc: "CSV file with high-level annotations of each high-confidence contig from cell-associated barcodes"

  loupe_vdj_browser_track:
    type: File?
    outputSource: aggregate_counts/loupe_vdj_browser_track
    label: "Loupe V(D)J Browser visualization and analysis file"
    doc: "Loupe V(D)J Browser visualization and analysis file"

  airr_rearrangement_tsv:
    type: File?
    outputSource: aggregate_counts/airr_rearrangement_tsv
    label: "Annotated contigs and consensus sequences of V(D)J rearrangements in the AIRR format"
    doc: |
      Annotated contigs and consensus sequences of V(D)J
      rearrangements in the AIRR format. It includes only
      viable cells identified by both V(D)J and RNA algorithms.

  compressed_html_data_folder:
    type: File
    outputSource: compress_html_data_folder/compressed_folder
    label: "Compressed folder with CellBrowser formatted results"
    doc: "Compressed folder with CellBrowser formatted results"

  html_data_folder:
    type: Directory
    outputSource: cellbrowser_build/html_data
    label: "Folder with not compressed CellBrowser formatted results"
    doc: "Folder with not compressed CellBrowser formatted results"

  cellbrowser_report:
    type: File
    outputSource: cellbrowser_build/index_html_file
    label: "CellBrowser formatted Cellranger report"
    doc: "CellBrowser formatted Cellranger report"
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  aggregate_counts_stdout_log:
    type: File
    outputSource: aggregate_counts/stdout_log
    label: "stdout log generated by cellranger aggr"
    doc: "stdout log generated by cellranger aggr"

  aggregate_counts_stderr_log:
    type: File
    outputSource: aggregate_counts/stderr_log
    label: "stderr log generated by cellranger aggr"
    doc: "stderr log generated by cellranger aggr"


steps:

  aggregate_counts:
    run: ../tools/cellranger-aggr.cwl
    in:
      molecule_info_h5: molecule_info_h5
      filtered_data_folder: filtered_data_folder
      gem_well_labels: gem_well_labels
      normalization_mode: normalization_mode
      clonotype_grouping: clonotype_grouping
      threads: threads
      memory_limit: memory_limit
      virt_memory_limit: memory_limit
    out:
    - web_summary_report
    - metrics_summary_report_json
    - secondary_analysis_report_folder
    - filtered_feature_bc_matrix_folder
    - filtered_feature_bc_matrix_h5
    - aggregation_metadata
    - grouping_data
    - loupe_browser_track
    - clonotypes_csv
    - consensus_sequences_fasta
    - consensus_annotations_csv
    - filtered_contig_annotations_csv
    - loupe_vdj_browser_track
    - airr_rearrangement_tsv
    - stdout_log
    - stderr_log

  compress_filtered_feature_bc_matrix_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: aggregate_counts/filtered_feature_bc_matrix_folder
    out:
    - compressed_folder

  compress_secondary_analysis_report_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: aggregate_counts/secondary_analysis_report_folder
    out:
    - compressed_folder

  cellbrowser_build:
    run: ../tools/cellbrowser-build-cellranger.cwl
    in:
      secondary_analysis_report_folder: aggregate_counts/secondary_analysis_report_folder
      filtered_feature_bc_matrix_folder: aggregate_counts/filtered_feature_bc_matrix_folder
      aggregation_metadata: aggregate_counts/aggregation_metadata
    out:
    - html_data
    - index_html_file

  compress_html_data_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: cellbrowser_build/html_data
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Cell Ranger Aggregate (RNA, RNA+VDJ)"
s:name: "Cell Ranger Aggregate (RNA, RNA+VDJ)"
s:alternateName: "Combines outputs from multiple runs of either Cell Ranger Count (RNA) or Cell Ranger Count (RNA+VDJ) pipelines"

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
  Cell Ranger Aggregate (RNA, RNA+VDJ)

  Combines outputs from multiple runs of either “Cell Ranger Count (RNA)”
  or “Cell Ranger Count (RNA+VDJ)” pipelines. The results of this workflow
  are primarily used in “Single-Cell RNA-Seq Filtering Analysis” and
  “Single-Cell Immune Profiling Analysis” pipelines.
