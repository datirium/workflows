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
    doc: "Cell Ranger ARC generated genome indices folder"
    'sd:upstreamSource': "genome_indices/arc_indices_folder"
    'sd:localLabel': true

  gex_fastq_file_r1:
    type:
    - File
    - type: array
      items: File
    format: "http://edamontology.org/format_1930"
    label: "GEX FASTQ file R1 (optionally compressed)"
    doc: "GEX FASTQ file R1 (optionally compressed)"

  gex_fastq_file_r2:
    type:
    - File
    - type: array
      items: File
    format: "http://edamontology.org/format_1930"
    label: "GEX FASTQ file R2 (optionally compressed)"
    doc: "GEX FASTQ file R2 (optionally compressed)"

  atac_fastq_file_r1:
    type:
    - File
    - type: array
      items: File
    format: "http://edamontology.org/format_1930"
    label: "ATAC FASTQ file R1 (optionally compressed)"
    doc: "ATAC FASTQ file R1 (optionally compressed)"

  atac_fastq_file_r2:
    type:
    - File
    - type: array
      items: File
    format: "http://edamontology.org/format_1930"
    label: "ATAC FASTQ file R2 (optionally compressed)"
    doc: "ATAC FASTQ file R2 (optionally compressed)"

  atac_fastq_file_r3:
    type:
    - File
    - type: array
      items: File
    format: "http://edamontology.org/format_1930"
    label: "ATAC FASTQ file R3 (optionally compressed)"
    doc: "ATAC FASTQ file R3 (optionally compressed)"

  exclude_introns:
    type: boolean?
    default: false
    label: "Disable counting of intronic reads"
    doc: |
      Disable counting of intronic reads. In this mode, only reads that are exonic
      and compatible with annotated splice junctions in the reference are counted.
      Note: using this mode will reduce the UMI counts in the feature-barcode matrix
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
    label: "Genome Type"
    doc: |
      Maximum memory used (GB).
      The same as was used for generating indices.
      The same will be applied to virtual memory
    'sd:upstreamSource': "genome_indices/memory_limit"
    'sd:localLabel': true


outputs:

  fastqc_report_gex_fastq_r1:
    type: File
    outputSource: run_fastqc_for_gex_fastq_r1/html_file
    label: "FastqQC report for GEX FASTQ file R1"
    doc: |
      FastqQC report for GEX FASTQ file R1
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  fastqc_report_gex_fastq_r2:
    type: File
    outputSource: run_fastqc_for_gex_fastq_r2/html_file
    label: "FastqQC report for GEX FASTQ file R2"
    doc: |
      FastqQC report for GEX FASTQ file R2
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  fastqc_report_atac_fastq_r1:
    type: File
    outputSource: run_fastqc_for_atac_fastq_r1/html_file
    label: "FastqQC report for ATAC FASTQ file R1"
    doc: |
      FastqQC report for ATAC FASTQ file R1
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  fastqc_report_atac_fastq_r2:
    type: File
    outputSource: run_fastqc_for_atac_fastq_r2/html_file
    label: "FastqQC report for ATAC FASTQ file R2"
    doc: |
      FastqQC report for ATAC FASTQ file R2
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  fastqc_report_atac_fastq_r3:
    type: File
    outputSource: run_fastqc_for_atac_fastq_r3/html_file
    label: "FastqQC report for ATAC FASTQ file R3"
    doc: |
      FastqQC report for ATAC FASTQ file R3
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

  barcode_metrics_report:
    type: File
    outputSource: generate_counts_matrix/barcode_metrics_report
    label: "ATAC and GEX barcode metrics in CSV format"
    doc: |
      ATAC and GEX read count summaries generated for every
      barcode observed in the experiment. The columns contain
      the paired ATAC and Gene Expression barcode sequences,
      ATAC and Gene Expression QC metrics for that barcode,
      as well as whether this barcode was identified as a
      cell-associated partition by the pipeline.

  gex_possorted_genome_bam_bai:
    type: File
    outputSource: generate_counts_matrix/gex_possorted_genome_bam_bai
    label: "Aligned to the genome indexed reads GEX BAM+BAI files"
    doc: |
      GEX position-sorted reads aligned to the genome and transcriptome
      annotated with barcode information in BAM format

  atac_possorted_genome_bam_bai:
    type: File
    outputSource: generate_counts_matrix/atac_possorted_genome_bam_bai
    label: "Aligned to the genome indexed reads ATAC BAM+BAI files"
    doc: |
      ATAC position-sorted reads aligned to the genome annotated with
      barcode information in BAM format

  filtered_feature_bc_matrix_folder:
    type: File
    outputSource: compress_filtered_feature_bc_matrix_folder/compressed_folder
    label: "Compressed folder with filtered feature-barcode matrices"
    doc: |
      Filtered feature barcode matrix stored as a CSC sparse matrix in MEX format.
      The rows consist of all the gene and peak features concatenated together
      (identical to raw feature barcode matrix) and the columns are restricted to
      those barcodes that are identified as cells.

  filtered_feature_bc_matrix_h5:
    type: File
    outputSource: generate_counts_matrix/filtered_feature_bc_matrix_h5
    label: "Filtered feature-barcode matrices in HDF5 format"
    doc: |
      Filtered feature barcode matrix stored as a CSC sparse matrix in hdf5 format.
      The rows consist of all the gene and peak features concatenated together
      (identical to raw feature barcode matrix) and the columns are restricted to
      those barcodes that are identified as cells.

  raw_feature_bc_matrices_folder:
    type: File
    outputSource: compress_raw_feature_bc_matrices_folder/compressed_folder
    label: "Compressed folder with unfiltered feature-barcode matrices"
    doc: |
      Raw feature barcode matrix stored as a CSC sparse matrix in MEX format.
      The rows consist of all the gene and peak features concatenated together
      and the columns consist of all observed barcodes with non-zero signal for
      either ATAC or gene expression.

  raw_feature_bc_matrices_h5:
    type: File
    outputSource: generate_counts_matrix/raw_feature_bc_matrices_h5
    label: "Unfiltered feature-barcode matrices in HDF5 format"
    doc: |
      Raw feature barcode matrix stored as a CSC sparse matrix in hdf5 format.
      The rows consist of all the gene and peak features concatenated together
      and the columns consist of all observed barcodes with non-zero signal for
      either ATAC or gene expression.

  secondary_analysis_report_folder:
    type: File
    outputSource: compress_secondary_analysis_report_folder/compressed_folder
    label: "Compressed folder with secondary analysis results"
    doc: |
      Various secondary analyses that utilize the ATAC data, the GEX data, and their
      linkage: dimensionality reduction and clustering results for the ATAC and GEX
      data, differential expression, and differential accessibility for all clustering
      results above and linkage between ATAC and GEX data.

  gex_molecule_info_h5:
    type: File
    outputSource: generate_counts_matrix/gex_molecule_info_h5
    label: "GEX molecule-level information for aggregating samples into larger datasets"
    doc: |
      Count and barcode information for every GEX molecule observed in the experiment
      in hdf5 format

  loupe_browser_track:
    type: File
    outputSource: generate_counts_matrix/loupe_browser_track
    label: "Loupe Browser visualization file with all the analysis outputs"
    doc: |
      Loupe Browser visualization file with all the analysis outputs

  atac_fragments_file:
    type: File
    outputSource: generate_counts_matrix/atac_fragments_file
    label: "Count and barcode information for every ATAC fragment in TSV format"
    doc: |
      Count and barcode information for every ATAC fragment observed in
      the experiment in TSV format.
  
  atac_peaks_bed_file:
    type: File
    outputSource: generate_counts_matrix/atac_peaks_bed_file
    label: "Identified peaks in BED format"
    doc: |
      Locations of open-chromatin regions identified in this sample.
      These regions are referred to as "peaks".

  atac_cut_sites_bigwig_file:
    type: File
    outputSource: generate_counts_matrix/atac_cut_sites_bigwig_file
    label: "Observed transposition sites in bigWig format"
    doc: |
      Genome track of observed transposition sites in the experiment
      smoothed at a resolution of 400 bases in BIGWIG format.
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "ATAC cut sites"
        height: 120

  atac_peak_annotation_file:
    type: File
    outputSource: generate_counts_matrix/atac_peak_annotation_file
    label: "Annotations of peaks based on genomic proximity in TSV format"
    doc: |
      Annotations of peaks based on genomic proximity alone.
      Note that these are not functional annotations and they
      do not make use of linkage with GEX data.

  generate_counts_matrix_stdout_log:
    type: File
    outputSource: generate_counts_matrix/stdout_log
    label: stdout log generated by cellranger-arc count
    doc: |
      stdout log generated by cellranger-arc count

  generate_counts_matrix_stderr_log:
    type: File
    outputSource: generate_counts_matrix/stderr_log
    label: stderr log generated by cellranger-arc count
    doc: |
      stderr log generated by cellranger-arc count

  collected_statistics:
    type: File
    outputSource: collect_statistics/collected_statistics
    label: "Collected statistics in Markdown format"
    doc: "Collected statistics in Markdown format"
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

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

  extract_gex_fastq_r1:
    run: ../tools/extract-fastq.cwl
    in:
      output_prefix:  
        default: "gex_read_1"
      compressed_file: gex_fastq_file_r1
    out:
    - fastq_file

  extract_gex_fastq_r2:
    run: ../tools/extract-fastq.cwl
    in:
      output_prefix:  
        default: "gex_read_2"
      compressed_file: gex_fastq_file_r2
    out:
    - fastq_file

  extract_atac_fastq_r1:
    run: ../tools/extract-fastq.cwl
    in:
      output_prefix: 
        default: "atac_read_1"
      compressed_file: atac_fastq_file_r1
    out:
    - fastq_file

  extract_atac_fastq_r2:
    run: ../tools/extract-fastq.cwl
    in:
      output_prefix:
        default: "atac_read_2"
      compressed_file: atac_fastq_file_r2
    out:
    - fastq_file

  extract_atac_fastq_r3:
    run: ../tools/extract-fastq.cwl
    in:
      output_prefix: 
        default: "atac_read_3"
      compressed_file: atac_fastq_file_r3
    out:
    - fastq_file


  run_fastqc_for_gex_fastq_r1:
    run: ../tools/fastqc.cwl
    in:
      reads_file: extract_gex_fastq_r1/fastq_file
      threads: threads
    out:
    - html_file

  run_fastqc_for_gex_fastq_r2:
    run: ../tools/fastqc.cwl
    in:
      reads_file: extract_gex_fastq_r2/fastq_file
      threads: threads
    out:
    - html_file

  run_fastqc_for_atac_fastq_r1:
    run: ../tools/fastqc.cwl
    in:
      reads_file: extract_atac_fastq_r1/fastq_file
      threads: threads
    out:
    - html_file

  run_fastqc_for_atac_fastq_r2:
    run: ../tools/fastqc.cwl
    in:
      reads_file: extract_atac_fastq_r2/fastq_file
      threads: threads
    out:
    - html_file

  run_fastqc_for_atac_fastq_r3:
    run: ../tools/fastqc.cwl
    in:
      reads_file: extract_atac_fastq_r3/fastq_file
      threads: threads
    out:
    - html_file


  generate_counts_matrix:
    run: ../tools/cellranger-arc-count.cwl
    in:
      gex_fastq_file_r1: extract_gex_fastq_r1/fastq_file
      gex_fastq_file_r2: extract_gex_fastq_r2/fastq_file
      atac_fastq_file_r1: extract_atac_fastq_r1/fastq_file
      atac_fastq_file_r2: extract_atac_fastq_r2/fastq_file
      atac_fastq_file_r3: extract_atac_fastq_r3/fastq_file
      indices_folder: indices_folder
      exclude_introns: exclude_introns
      threads: threads
      memory_limit: memory_limit
      virt_memory_limit: memory_limit
    out:
    - web_summary_report
    - metrics_summary_report
    - barcode_metrics_report
    - gex_possorted_genome_bam_bai
    - atac_possorted_genome_bam_bai
    - filtered_feature_bc_matrix_folder
    - filtered_feature_bc_matrix_h5
    - raw_feature_bc_matrices_folder
    - raw_feature_bc_matrices_h5
    - secondary_analysis_report_folder
    - gex_molecule_info_h5
    - loupe_browser_track
    - atac_fragments_file
    - atac_peaks_bed_file
    - atac_cut_sites_bigwig_file
    - atac_peak_annotation_file
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
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: rackspacedot/python37
      inputs:
        script:
          type: string?
          default: |
            #!/usr/bin/env python3
            import sys, csv
            with open(sys.argv[1], "r") as input_stream:
              with open("collected_statistics.md", "w") as output_stream:
                output_stream.write("### Cell Ranger ARC Statistics\n")
                keys, values = None, None
                for i, row in enumerate(csv.reader(input_stream)):
                  if i==0:
                    keys = row
                  else:
                    values = row
                for k,v in zip(keys, values):
                  output_stream.write("- "+k+": "+v+"\n")
          inputBinding:
            position: 5
        metrics_summary_report:
          type: File
          inputBinding:
            position: 6
      outputs:
        collected_statistics:
          type: File
          outputBinding:
            glob: "*"
      baseCommand: ["python3", "-c"]
    in:
      metrics_summary_report: generate_counts_matrix/metrics_summary_report
    out:
    - collected_statistics

  cellbrowser_build:
    run: ../tools/cellbrowser-build-cellranger-arc.cwl
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

s:name: "Cell Ranger ARC Count Gene Expression + ATAC"
label: "Cell Ranger ARC Count Gene Expression + ATAC"
s:alternateName: "Counts ATAC and gene expression reads from a single 10x Genomics Cell Ranger Multiome ATAC + Gene Expression library"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/cellranger-arc-count.cwl
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
  Cell Ranger ARC Count Gene Expression + ATAC
  ============================================