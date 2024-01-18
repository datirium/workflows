cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  genome_indices:
    - "cellranger-mkref.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  indices_folder:
    type: Directory
    label: "Genome Type"
    doc: "Cell Ranger generated genome indices folder"
    'sd:upstreamSource': "genome_indices/indices_folder"
    'sd:localLabel': true

  fastq_file_r1:
    type:
    - File
    - type: array
      items: File
    label: "FASTQ file(s) R1 (optionally compressed)"
    doc: "FASTQ file(s) R1 (optionally compressed)"

  fastq_file_r2:
    type:
    - File
    - type: array
      items: File
    label: "FASTQ file(s) R2 (optionally compressed)"
    doc: "FASTQ file(s) R2 (optionally compressed)"

  r1_length:
    type: int?
    default: null
    label: "Limit the length of the input R1 sequence"
    doc: "Limit the length of the input R1 sequence"
    'sd:layout':
      advanced: true

  r2_length:
    type: int?
    default: null
    label: "Limit the length of the input R2 sequence"
    doc: "Limit the length of the input R2 sequence"
    'sd:layout':
      advanced: true

  expect_cells:
    type: int?
    default: null
    label: "Expected number of recovered cells. If not provided - use auto-estimated"
    doc: "Expected number of recovered cells. If not provided - use auto-estimated"
    'sd:layout':
      advanced: true

  force_cells:
    type: int?
    default: null
    label: "Force pipeline to use this number of cells, bypassing the cell detection algorithm"
    doc: "Force pipeline to use this number of cells, bypassing the cell detection algorithm"
    'sd:layout':
      advanced: true

  exclude_introns:
    type: boolean?
    default: false
    label: "Do not count intronic reads for whole transcriptome gene expression data"
    doc: "Do not count intronic reads for whole transcriptome gene expression data"
    'sd:layout':
      advanced: true

  no_bam:
    type: boolean?
    default: true
    label: "Do not generate the BAM file"
    doc: "Do not generate the BAM file"
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
    default: 30
    label: "Genome Type"
    doc: |
      Maximum memory used (GB).
      The same as was used for generating indices.
      The same will be applied to virtual memory
    'sd:upstreamSource': "genome_indices/memory_limit"
    'sd:localLabel': true


outputs:

  fastqc_report_fastq_r1:
    type: File
    outputSource: run_fastqc_for_fastq_r1/html_file
    label: "FastqQC report for FASTQ file R1"
    doc: |
      FastqQC report for FASTQ file R1
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  fastqc_report_fastq_r2:
    type: File
    outputSource: run_fastqc_for_fastq_r2/html_file
    label: "FastqQC report for FASTQ file R2"
    doc: |
      FastqQC report for FASTQ file R2
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  web_summary_report:
    type: File
    outputSource: generate_counts_matrix/web_summary_report
    label: "Cell Ranger summary"
    doc: |
      Cell Ranger summary
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  metrics_summary_report:
    type: File
    outputSource: generate_counts_matrix/metrics_summary_report
    label: "Run summary metrics in CSV format"
    doc: |
      Run summary metrics in CSV format

  possorted_genome_bam_bai:
    type: File?
    outputSource: generate_counts_matrix/possorted_genome_bam_bai
    label: "Aligned to the genome indexed reads BAM+BAI files"
    doc: |
      Indexed reads aligned to the genome and transcriptome annotated
      with barcode information

  filtered_feature_bc_matrix_folder:
    type: File
    outputSource: compress_filtered_feature_bc_matrix_folder/compressed_folder
    label: "Compressed folder with filtered feature-barcode matrices"
    doc: |
      Compressed folder with filtered feature-barcode matrices containing only cellular barcodes in MEX format.
      When implemented, in Targeted Gene Expression samples, the non-targeted genes won't be present.

  filtered_feature_bc_matrix_h5:
    type: File
    outputSource: generate_counts_matrix/filtered_feature_bc_matrix_h5
    label: "Filtered feature-barcode matrices in HDF5 format"
    doc: |
      Filtered feature-barcode matrices containing only cellular barcodes in HDF5 format.
      When implemented, in Targeted Gene Expression samples, the non-targeted genes won't
      be present.
  
  raw_feature_bc_matrices_folder:
    type: File
    outputSource: compress_raw_feature_bc_matrices_folder/compressed_folder
    label: "Compressed folder with unfiltered feature-barcode matrices"
    doc: |
      Compressed folder with unfiltered feature-barcode matrices containing all barcodes in MEX format

  raw_feature_bc_matrices_h5:
    type: File
    outputSource: generate_counts_matrix/raw_feature_bc_matrices_h5
    label: "Unfiltered feature-barcode matrices in HDF5 format"
    doc: |
      Unfiltered feature-barcode matrices containing all barcodes in HDF5 format

  secondary_analysis_report_folder:
    type: File
    outputSource: compress_secondary_analysis_report_folder/compressed_folder
    label: "Compressed folder with secondary analysis results"
    doc: |
      Compressed folder with secondary analysis results including dimensionality reduction,
      cell clustering, and differential expression

  molecule_info_h5:
    type: File
    outputSource: generate_counts_matrix/molecule_info_h5
    label: "Molecule-level information for aggregating samples into larger datasets"
    doc: |
      Molecule-level information used by cellranger aggr to aggregate samples into
      larger datasets

  loupe_browser_track:
    outputSource: generate_counts_matrix/loupe_browser_track
    label: "Loupe Browser visualization and analysis file"
    type: File
    doc: |
      Loupe Browser visualization and analysis file

  collected_statistics_yaml:
    type: File
    outputSource: collect_statistics/collected_statistics_yaml
    label: "Collected statistics in YAML format"
    doc: "Collected statistics in YAML format"

  collected_statistics_md:
    type: File
    outputSource: collect_statistics/collected_statistics_md
    label: "Collected statistics in Markdown format"
    doc: "Collected statistics in Markdown format"
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  collected_statistics_tsv:
    type: File
    outputSource: collect_statistics/collected_statistics_tsv
    label: "Collected statistics in TSV format"
    doc: "Collected statistics in TSV format"
    'sd:visualPlugins':
    - tableView:
        vertical: true
        tab: 'Overview'

  generate_counts_matrix_stdout_log:
    type: File
    outputSource: generate_counts_matrix/stdout_log
    label: stdout log generated by cellranger count
    doc: |
      stdout log generated by cellranger count

  generate_counts_matrix_stderr_log:
    type: File
    outputSource: generate_counts_matrix/stderr_log
    label: stderr log generated by cellranger count
    doc: |
      stderr log generated by cellranger count

  compressed_html_data_folder:
    type: File
    outputSource: compress_html_data_folder/compressed_folder
    label: "Compressed folder with CellBrowser formatted results"
    doc: |
      Compressed folder with CellBrowser formatted results

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
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"


steps:

  extract_fastq_r1:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_r1
      output_prefix:
        default: "read_1"
    out:
    - fastq_file

  extract_fastq_r2:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_r2
      output_prefix:
        default: "read_2"
    out:
    - fastq_file

  run_fastqc_for_fastq_r1:
    run: ../tools/fastqc.cwl
    in:
      reads_file: extract_fastq_r1/fastq_file
      threads: threads
    out:
    - html_file

  run_fastqc_for_fastq_r2:
    run: ../tools/fastqc.cwl
    in:
      reads_file: extract_fastq_r2/fastq_file
      threads: threads
    out:
    - html_file

  generate_counts_matrix:
    run: ../tools/cellranger-count.cwl
    in:
      fastq_file_r1: extract_fastq_r1/fastq_file
      fastq_file_r2: extract_fastq_r2/fastq_file
      indices_folder: indices_folder
      r1_length: r1_length
      r2_length: r2_length
      expect_cells: expect_cells
      force_cells: force_cells
      no_bam: no_bam
      exclude_introns: exclude_introns
      threads: threads
      memory_limit: memory_limit
      virt_memory_limit: memory_limit
    out:
    - web_summary_report
    - metrics_summary_report
    - possorted_genome_bam_bai
    - filtered_feature_bc_matrix_folder
    - filtered_feature_bc_matrix_h5
    - raw_feature_bc_matrices_folder
    - raw_feature_bc_matrices_h5
    - secondary_analysis_report_folder
    - molecule_info_h5
    - loupe_browser_track
    - stdout_log
    - stderr_log

  compress_filtered_feature_bc_matrix_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: generate_counts_matrix/filtered_feature_bc_matrix_folder
    out:
    - compressed_folder

  compress_raw_feature_bc_matrices_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: generate_counts_matrix/raw_feature_bc_matrices_folder
    out:
    - compressed_folder

  compress_secondary_analysis_report_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: generate_counts_matrix/secondary_analysis_report_folder
    out:
    - compressed_folder

  collect_statistics:
    run: ../tools/collect-stats-sc-count.cwl
    in:
      metrics_summary_report: generate_counts_matrix/metrics_summary_report
    out:
    - collected_statistics_yaml
    - collected_statistics_tsv
    - collected_statistics_md

  cellbrowser_build:
    run: ../tools/cellbrowser-build-cellranger.cwl
    in:
      secondary_analysis_report_folder: generate_counts_matrix/secondary_analysis_report_folder
      filtered_feature_bc_matrix_folder: generate_counts_matrix/filtered_feature_bc_matrix_folder
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

s:name: "Cell Ranger Count (RNA)"
label: "Cell Ranger Count (RNA)"
s:alternateName: "Quantifies single-cell gene expression of the sequencing data from a single 10x Genomics library"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/single-cell-preprocess-cellranger.cwl
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
  Cell Ranger Count (RNA)

  Quantifies single-cell gene expression of the sequencing data
  from a single 10x Genomics library. The results of this workflow
  are primarily used in either “Single-Cell RNA-Seq Filtering
  Analysis” or “Cell Ranger Aggregate (RNA, RNA+VDJ)” pipelines.