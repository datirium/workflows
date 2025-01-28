cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


"sd:upstream":
  genome_indices:
  - "cellranger-mkref.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

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

  multiome_arc:
    type: boolean?
    default: false
    label: "scATAC-Seq files come from scMultiome experiment"
    doc: |
      Changes chemistry type parameter to indicate
      that scATAC-Seq data is part of the scMultiome
      experiment.

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

  fastq_file_r3:
    type:
    - File
    - type: array
      items: File
    label: "FASTQ file(s) R3 (optionally compressed)"
    doc: "FASTQ file(s) R3 (optionally compressed)"

  force_cells:
    type: int?
    default: null
    label: "Define the top N barcodes with the most ATAC fragments overlapping peaks as cells"
    doc: |
      Define the top N barcodes with the most ATAC fragments overlapping
      peaks as cells. N must be a positive integer <= 20,000. Please
      consult the documentation before using this option
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
    outputSource: generate_counts_matrix/web_summary_report
    label: "Cell Ranger summary"
    doc: |
      Run summary metrics and charts in HTML format
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  metrics_summary_report_json:
    type: File
    outputSource: generate_counts_matrix/metrics_summary_report_json
    label: "Run summary metrics in JSON format"
    doc: |
      Run summary metrics in JSON format

  metrics_summary_report_csv:
    type: File
    outputSource: generate_counts_matrix/metrics_summary_report_csv
    label: "Run summary metrics in CSV format"
    doc: |
      Run summary metrics in CSV format

  barcode_metrics_report:
    type: File
    outputSource: generate_counts_matrix/barcode_metrics_report
    label: "Per-barcode fragment counts & metrics in CSV format"
    doc: |
      Per-barcode fragment counts & metrics in CSV format

  possorted_genome_bam_bai:
    type: File
    outputSource: generate_counts_matrix/possorted_genome_bam_bai
    label: "ATAC reads"
    doc: |
      Genome track of ATAC reads aligned to
      the reference genome. Each read has
      a 10x Chromium cellular (associated
      with a 10x Genomics gel bead) barcode
      and mapping information stored in TAG
      fields.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "IGV Genome Browser"
        id: "igvbrowser"
        type: "alignment"
        format: "bam"
        name: "ATAC reads"
        displayMode: "SQUISHED"

  atac_fragments_file:
    type: File
    outputSource: generate_counts_matrix/fragments_file
    label: "Count and barcode information for every ATAC fragment in TSV format"
    doc: |
      Count and barcode information for every ATAC fragment observed
      in the experiment in TSV format

  peaks_bed_file:
    type: File
    outputSource: generate_counts_matrix/peaks_bed_file
    label: "ATAC peaks"
    doc: |
      Genome track of open-chromatin
      regions identified as peaks.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "IGV Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        name: "ATAC peaks"
        displayMode: "COLLAPSE"
        height: 40

  peak_annotation_file:
    type: File
    outputSource: generate_counts_matrix/peak_annotation_file
    label: "Annotations of peaks based on genomic proximity alone"
    doc: |
      Annotations of peaks based on genomic proximity alone

  cut_sites_bigwig_file:
    type: File
    outputSource: generate_counts_matrix/cut_sites_bigwig_file
    label: "ATAC transposition counts"
    doc: |
      Genome track of observed transposition
      sites in the experiment smoothed at a
      resolution of 400 bases.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "IGV Genome Browser"
        id: "igvbrowser"
        type: "wig"
        name: "ATAC transposition counts"
        height: 120

  # peak_motif_mapping_bed:
  #   type: File
  #   outputSource: generate_counts_matrix/peak_motif_mapping_bed
  #   label: "File with peak-motif associations in BED format"
  #   doc: |
  #     File with peak-motif associations in BED format

  filtered_feature_bc_matrix_folder:
    type: File
    outputSource: compress_filtered_feature_bc_matrix_folder/compressed_folder
    label: "Compressed folder with filtered peak-barcode matrices"
    doc: |
      Folder with filtered peak-barcode matrices containing only
      cellular barcodes in MEX format.

  filtered_feature_bc_matrix_h5:
    type: File
    outputSource: generate_counts_matrix/filtered_feature_bc_matrix_h5
    label: "Filtered peak-barcode matrices in HDF5 format"
    doc: |
      Filtered peak-barcode matrices containing only cellular
      barcodes in HDF5 format.

  # filtered_tf_bc_matrix_folder:
  #   type: File
  #   outputSource: compress_filtered_tf_bc_matrix_folder/compressed_folder
  #   label: "Compressed folder with filtered tf-barcode matrices"
  #   doc: |
  #     Folder with filtered tf-barcode matrices containing only cellular
  #     barcodes in MEX format.

  # filtered_tf_bc_matrix_h5:
  #   type: File
  #   outputSource: generate_counts_matrix/filtered_tf_bc_matrix_h5
  #   label: "Filtered tf-barcode matrices in HDF5 format"
  #   doc: |
  #     Filtered tf-barcode matrices containing only cellular
  #     barcodes in HDF5 format.

  raw_feature_bc_matrices_folder:
    type: File
    outputSource: compress_raw_feature_bc_matrices_folder/compressed_folder
    label: "Compressed folder with unfiltered peak-barcode matrices"
    doc: |
      Folder with unfiltered peak-barcode matrices containing
      all barcodes in MEX format

  raw_feature_bc_matrices_h5:
    type: File
    outputSource: generate_counts_matrix/raw_feature_bc_matrices_h5
    label: "Unfiltered peak-barcode matrices in HDF5 format"
    doc: |
      Unfiltered peak-barcode matrices containing all barcodes
      in HDF5 format

  secondary_analysis_report_folder:
    type: File
    outputSource: compress_secondary_analysis_report_folder/compressed_folder
    label: "Compressed folder with secondary analysis results"
    doc: |
      Folder with secondary analysis results

  loupe_browser_track:
    type: File
    outputSource: generate_counts_matrix/loupe_browser_track
    label: "Loupe Browser visualization and analysis file"
    doc: |
      Loupe Browser visualization and analysis file

  generate_counts_matrix_stdout_log:
    type: File
    outputSource: generate_counts_matrix/stdout_log
    label: stdout log generated by cellranger-atac count
    doc: |
      stdout log generated by cellranger-atac count

  generate_counts_matrix_stderr_log:
    type: File
    outputSource: generate_counts_matrix/stderr_log
    label: stderr log generated by cellranger-atac count
    doc: |
      stderr log generated by cellranger-atac count

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
    "sd:visualPlugins":
    - markdownView:
        tab: "Overview"

  collected_statistics_tsv:
    type: File
    outputSource: collect_statistics/collected_statistics_tsv
    label: "Collected statistics in TSV format"
    doc: "Collected statistics in TSV format"
    "sd:visualPlugins":
    - tableView:
        vertical: true
        tab: "Overview"

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

  extract_fastq_r3:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_r3
      output_prefix:
        default: "read_3"
    out:
    - fastq_file

  generate_counts_matrix:
    run: ../tools/cellranger-atac-count.cwl
    in:
      fastq_file_r1: extract_fastq_r1/fastq_file
      fastq_file_r2: extract_fastq_r2/fastq_file
      fastq_file_r3: extract_fastq_r3/fastq_file
      indices_folder: indices_folder
      force_cells: force_cells
      chemistry:
        source: multiome_arc
        valueFrom: $(self?"ARC-v1":"null")
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
    - possorted_genome_bam_bai
    - fragments_file
    - peaks_bed_file
    - peak_annotation_file
    - cut_sites_bigwig_file
    # - peak_motif_mapping_bed
    - filtered_feature_bc_matrix_folder
    - filtered_feature_bc_matrix_h5
    # - filtered_tf_bc_matrix_folder
    # - filtered_tf_bc_matrix_h5
    - raw_feature_bc_matrices_folder
    - raw_feature_bc_matrices_h5
    - secondary_analysis_report_folder
    - loupe_browser_track
    - stdout_log
    - stderr_log

  compress_filtered_feature_bc_matrix_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: generate_counts_matrix/filtered_feature_bc_matrix_folder
    out:
    - compressed_folder

  # compress_filtered_tf_bc_matrix_folder:
  #   run: ../tools/tar-compress.cwl
  #   in:
  #     folder_to_compress: generate_counts_matrix/filtered_tf_bc_matrix_folder
  #   out:
  #   - compressed_folder

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
    run: ../tools/collect-stats-sc-atac-count.cwl
    in:
      metrics_summary_report: generate_counts_matrix/metrics_summary_report_csv
    out:
    - collected_statistics_yaml
    - collected_statistics_tsv
    - collected_statistics_md

  cellbrowser_build:
    run: ../tools/cellbrowser-build-cellranger-atac.cwl
    in:
      secondary_analysis_report_folder: generate_counts_matrix/secondary_analysis_report_folder
      filtered_feature_bc_matrix_folder: generate_counts_matrix/filtered_feature_bc_matrix_folder
      annotation_gtf_file: annotation_gtf_file
    out:
    - html_data
    - index_html_file


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Cell Ranger Count (ATAC)"
s:name: "Cell Ranger Count (ATAC)"
s:alternateName: "Cell Ranger Count (ATAC)"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/cellranger-atac-count.cwl
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
  Cell Ranger Count (ATAC)

  Quantifies single-cell chromatin accessibility of the
  sequencing data from a single 10x Genomics library.
  The results of this workflow are used in either the
  “Single-Cell ATAC-Seq Filtering Analysis” or “Cell
  Ranger Aggregate (ATAC)” pipeline.