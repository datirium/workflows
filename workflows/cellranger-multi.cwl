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
  vdj_indices:
  - "cellranger-mkvdjref.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  rna_indices_folder:
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
    "sd:upstreamSource": "genome_indices/indices_folder"
    "sd:localLabel": true

  memory_limit:
    type: int?
    default: 20
    "sd:upstreamSource": "genome_indices/memory_limit"

  vdj_indices_folder:
    type: Directory
    label: "Cell Ranger Reference VDJ Sample"
    doc: |
      Any "Cell Ranger Reference VDJ Sample"
      that builds a reference genome of a
      selected species for V(D)J contigs
      assembly and clonotype calling. This
      sample can be obtained from "Cell
      Ranger Reference (VDJ)" pipeline.
    "sd:upstreamSource": "vdj_indices/indices_folder"
    "sd:localLabel": true

  rna_fastq_file_r1:
    type:
    - File
    - type: array
      items: File
    format: "http://edamontology.org/format_1930"
    label: "RNA FASTQ file(s) R1 (optionally compressed)"
    doc: "RNA FASTQ file(s) R1 (optionally compressed)"

  rna_fastq_file_r2:
    type:
    - File
    - type: array
      items: File
    format: "http://edamontology.org/format_1930"
    label: "RNA FASTQ file(s) R2 (optionally compressed)"
    doc: "RNA FASTQ file(s) R2 (optionally compressed)"

  vdj_fastq_file_r1:
    type:
    - File
    - type: array
      items: File
    format: "http://edamontology.org/format_1930"
    label: "V(D)J FASTQ file(s) R1 (optionally compressed)"
    doc: "V(D)J FASTQ file(s) R1 (optionally compressed)"

  vdj_fastq_file_r2:
    type:
    - File
    - type: array
      items: File
    format: "http://edamontology.org/format_1930"
    label: "V(D)J FASTQ file(s) R2 (optionally compressed)"
    doc: "V(D)J FASTQ file(s) R2 (optionally compressed)"

  vdj_chain_type:
    type:
    - "null"
    - type: enum
      name: "chain_type"
      symbols:
      - "VDJ"
      - "VDJ-T"
      - "VDJ-B"
      - "VDJ-T-GD"
    default: "VDJ"
    label: "V(D)J chain type. Use VDJ for auto detection."
    doc: |
      V(D)J chain type. Setting to VDJ will auto-detect the chain type.
      Auto-detection does not work for TRG/D (gamma-delta) chains.
      Note that gamma-delta analysis is enabled but the algorithm has
      not been tested extensively.
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
    outputSource: cellranger_multi/web_summary_report
    label: "Gene Expression and V(D)J Repertoire Profiling"
    doc: |
      Run summary metrics and charts in HTML format
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  metrics_summary_report:
    type: File
    outputSource: convert_metrics_summary_report_to_tsv/output_file
    label: "Run summary metrics in TSV format"
    doc: |
      Run summary metrics in TSV format
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "QC metrics"
        Title: "QC metrics"

  possorted_genome_bam_bai:
    type: File
    outputSource: cellranger_multi/possorted_genome_bam_bai
    label: "Unaligned and aligned to the genome and transcriptome indexed reads"
    doc: |
      Indexed RNA BAM file containing position-sorted reads aligned
      to the genome and transcriptome, as well as unaligned reads.

  filtered_feature_bc_matrix_folder:
    type: File
    outputSource: compress_filtered_feature_bc_matrix_folder/compressed_folder
    label: "Filtered feature-barcode matrices in MEX format"
    doc: |
      Folder with filtered feature-barcode matrices containing only cellular
      barcodes in MEX format. Each element of the matrix is the number of UMIs
      associated with a feature (row) and a barcode (column).

  filtered_feature_bc_matrix_h5:
    type: File
    outputSource: cellranger_multi/filtered_feature_bc_matrix_h5
    label: "Filtered feature-barcode matrices in HDF5 format"
    doc: |
      Filtered feature-barcode matrices containing only cellular
      barcodes in HDF5 format. Each element of the matrix is the
      number of UMIs associated with a feature (row) and a
      barcode (column).

  raw_feature_bc_matrices_folder:
    type: File
    outputSource: compress_raw_feature_bc_matrices_folder/compressed_folder
    label: "Unfiltered feature-barcode matrices in MEX format"
    doc: |
      Folder with unfiltered feature-barcode matrices containing all barcodes
      in MEX format. Each element of the matrix is the number of UMIs associated
      with a feature (row) and a barcode (column).

  raw_feature_bc_matrices_h5:
    type: File
    outputSource: cellranger_multi/raw_feature_bc_matrices_h5
    label: "Unfiltered feature-barcode matrices in HDF5 format"
    doc: |
      Unfiltered feature-barcode matrices containing all barcodes in HDF5 format.
      Each element of the matrix is the number of UMIs associated with a feature
      (row) and a barcode (column).

  secondary_analysis_report_folder:
    type: File
    outputSource: compress_secondary_analysis_report_folder/compressed_folder
    label: "Folder with secondary analysis of RNA data"
    doc: |
      Folder with secondary analysis of RNA data including dimensionality
      reduction, cell clustering, and differential expression

  loupe_browser_track:
    type: File
    outputSource: cellranger_multi/loupe_browser_track
    label: "Loupe Browser visualization and analysis file"
    doc: |
      Loupe Browser visualization and analysis file

  all_contig_reads_bam_bai:
    type: File
    outputSource: cellranger_multi/all_contig_reads_bam_bai
    label: "Indexed V(D)J BAM file with reads aligned to ALL assembled contigs, per cell barcode"
    doc: |
      Indexed V(D)J BAM file with reads aligned to ALL assembled contigs, per cell barcode.
      This file demonstrates how the reads and UMIs support the assembled contigs within
      a cell barcode. Reads are not aligned across cell barcode boundaries. Please note
      that this BAM excludes reads whose barcodes don't match the whitelist, so it is not
      suitable as an archive of every single input read.
      This file includes reads from all cells barcodes identified by V(D)J algorithm including
      those ones that will be later discarded as non-viable cells by V(D)J algorithm and those
      barcodes that will be later removed after overlapping with cells called by RNA algorithm.

  all_contig_sequences_fasta:
    type: File
    outputSource: cellranger_multi/all_contig_sequences_fasta
    label: "FASTA format sequence for ALL assembled contigs in the V(D)J library"
    doc: |
      FASTA format sequence for ALL assembled contigs in the V(D)J library.
      This file includes both productive and non-productive contigs with high and low confidence
      assembled for all identified cells barcodes including those ones that will be later discarded
      as non-viable cells by V(D)J algorithm or after overlapping with cells called by RNA algorithm.

  all_contig_annotations_bed:
    type: File
    outputSource: cellranger_multi/all_contig_annotations_bed
    label: "BED file with high-level and detailed annotations of ALL assembled contigs (from cell and background barcodes)"
    doc: |
      BED file with high-level and detailed annotations of ALL assembled contigs (from cell and
      background barcodes). Used for further investigation into why some contigs were filtered
      out. This file includes both productive and non-productive contigs with high and low
      confidence assembled for all identified cells barcodes including those ones that will be
      later discarded as non-viable cells by V(D)J algorithm or after overlapping with cells
      called by RNA algorithm.

  all_contig_annotations_csv:
    type: File
    outputSource: cellranger_multi/all_contig_annotations_csv
    label: "CSV file with high-level and detailed annotations of ALL assembled contigs (from cell and background barcodes)"
    doc: |
      CSV file with high-level and detailed annotations of ALL assembled contigs (from cell and
      background barcodes). Used for further investigation into why some contigs were filtered
      out. This file includes both productive and non-productive contigs with high and low
      confidence assembled for all identified cells barcodes including those ones that will be
      later discarded as non-viable cells by V(D)J algorithm or after overlapping with cells
      called by RNA algorithm.

  airr_rearrangement_tsv:
    type: File
    outputSource: cellranger_multi/airr_rearrangement_tsv
    label: "Annotated contigs and consensus sequences of V(D)J rearrangements in the AIRR format"
    doc: |
      Annotated contigs and consensus sequences of V(D)J rearrangements
      in the AIRR format. It includes only viable cells identified by
      both V(D)J and RNA algorithms.

  clonotypes_tsv:
    type: File
    outputSource: convert_clonotypes_csv_to_tsv/output_file
    label: "TSV file with high-level descriptions of each clonotype"
    doc: |
      TSV file with high-level descriptions of each clonotype. During the clonotype
      grouping stage, cell barcodes are placed in groups called clonotypes. Only viable
      cells identified by both V(D)J and RNA algorithms are used. Each clonotype consists
      of all descendants of a single, fully rearranged common ancestor, as approximated
      computationally. During this process, some cell barcodes are flagged as likely
      artifacts and filtered out, meaning that they are no longer called as cells.
      However, as clonotype grouping stage is hapenning before forming the final version
      of files in the per_sample_outs folder, the reported cells number won't be affected.
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "V(D)J clonotypes"
        Title: "V(D)J clonotypes"

  germline_contigs_bam_bai:
    type: File
    outputSource: cellranger_multi/germline_contigs_bam_bai
    label: "Indexed V(D)J BAM file with contigs aligned to concatenated germline segments"
    doc: |
      Indexed V(D)J BAM file with contigs aligned to concatenated germline
      segments. For each clonotype consensus, the reference sequence is the
      annotated germline segments concatenated together. This file shows how
      both the per-cell contigs and the clonotype consensus contig relate to
      the germline reference. Useful for revealing polymorphisms, somatic
      mutations, and recombination-induced differences such as non-templated
      nucleotide additions.

  germline_sequences_fasta:
    type: File
    outputSource: cellranger_multi/germline_sequences_fasta
    label: "Concatenated V(D)J reference segments for the segments detected on each consensus sequence"
    doc: |
      Concatenated V(D)J reference segments for the segments detected on each
      consensus sequence. These serve as an approximate reference for each
      consensus sequence.

  consensus_contigs_bam_bai:
    type: File
    outputSource: cellranger_multi/consensus_contigs_bam_bai
    label: "Indexed V(D)J BAM file with contigs aligned to clonotype consensus"
    doc: |
      Indexed V(D)J BAM file with contigs aligned to clonotype consensus.
      Each "reference" sequence is a clonotype consensus sequence, and each
      record is an alignment of a single cell's contig against this consensus.
      This file shows, for a clonotype consensus sequences, how the constituent
      per-cell assemblies support the consensus.

  consensus_sequences_fasta:
    type: File
    outputSource: cellranger_multi/consensus_sequences_fasta
    label: "The consensus sequence of each assembled contig"
    doc: |
      The consensus sequence of each assembled contig.

  consensus_annotations_csv:
    type: File
    outputSource: cellranger_multi/consensus_annotations_csv
    label: "CSV file with high-level and detailed annotations of each clonotype consensus sequence"
    doc: |
      CSV file with high-level and detailed annotations of each clonotype
      consensus sequence.

  filtered_contig_annotations_csv:
    type: File
    outputSource: cellranger_multi/filtered_contig_annotations_csv
    label: "CSV file with high-level annotations of each high-confidence contig from cell-associated barcodes"
    doc: |
      CSV file with high-level annotations of each high-confidence contig from
      cell-associated barcodes. This is a subset of all_contig_annotations.csv.

  filtered_contig_sequences_fasta:
    type: File
    outputSource: cellranger_multi/filtered_contig_sequences_fasta
    label: "FASTA format sequence for only high-confidence contigs in cell barcodes"
    doc: |
      FASTA format sequence for only high-confidence contigs in cell barcodes.

  loupe_vdj_browser_track:
    type: File
    outputSource: cellranger_multi/loupe_vdj_browser_track
    label: "Loupe V(D)J Browser visualization and analysis file"
    doc: |
      Loupe V(D)J Browser visualization and analysis file

  filtered_data_folder:
    type: Directory
    outputSource: cellranger_multi/filtered_data_folder
    label: "Folder containing filtered data, i.e., only cell-associated barcodes"
    doc: |
      Folder containing filtered data, i.e., only cell-associated barcodes.
      Used by cellranger aggr to aggregate samples for joint analysis.

  html_data_folder:
    type: Directory
    outputSource: cellbrowser_build/html_data
    label: "Folder with not compressed CellBrowser formatted results"
    doc: |
      Folder with not compressed CellBrowser formatted results

  cellbrowser_report:
    type: File
    outputSource: cellbrowser_build/index_html_file
    label: "UCSC Cell Browser"
    doc: |
      CellBrowser formatted Cellranger report
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  cellranger_multi_stdout_log:
    type: File
    outputSource: cellranger_multi/stdout_log
    label: stdout log generated by cellranger multi
    doc: |
      stdout log generated by cellranger multi

  cellranger_multi_stderr_log:
    type: File
    outputSource: cellranger_multi/stderr_log
    label: stderr log generated by cellranger multi
    doc: |
      stderr log generated by cellranger multi


steps:

  extract_rna_fastq_r1:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: rna_fastq_file_r1
      output_prefix:
        default: "rna_read_1"
    out:
    - fastq_file

  extract_rna_fastq_r2:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: rna_fastq_file_r2
      output_prefix:
        default: "rna_read_2"
    out:
    - fastq_file

  extract_vdj_fastq_r1:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: vdj_fastq_file_r1
      output_prefix:
        default: "vdj_read_1"
    out:
    - fastq_file

  extract_vdj_fastq_r2:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: vdj_fastq_file_r2
      output_prefix:
        default: "vdj_read_2"
    out:
    - fastq_file

  cellranger_multi:
    run: ../tools/cellranger-multi.cwl
    in:
      rna_fastq_file_r1: extract_rna_fastq_r1/fastq_file
      rna_fastq_file_r2: extract_rna_fastq_r2/fastq_file
      vdj_fastq_file_r1: extract_vdj_fastq_r1/fastq_file
      vdj_fastq_file_r2: extract_vdj_fastq_r2/fastq_file
      rna_indices_folder: rna_indices_folder
      vdj_indices_folder: vdj_indices_folder
      vdj_chain_type: vdj_chain_type
      threads:
        source: threads
        valueFrom: $(parseInt(self))
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
    - loupe_browser_track
    - all_contig_reads_bam_bai
    - all_contig_sequences_fasta
    - all_contig_annotations_bed
    - all_contig_annotations_csv
    - airr_rearrangement_tsv
    - clonotypes_csv
    - germline_contigs_bam_bai
    - germline_sequences_fasta
    - consensus_contigs_bam_bai
    - consensus_sequences_fasta
    - consensus_annotations_csv
    - filtered_contig_annotations_csv
    - filtered_contig_sequences_fasta
    - loupe_vdj_browser_track
    - filtered_data_folder
    - stdout_log
    - stderr_log

  compress_filtered_feature_bc_matrix_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: cellranger_multi/filtered_feature_bc_matrix_folder
    out:
    - compressed_folder

  compress_raw_feature_bc_matrices_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: cellranger_multi/raw_feature_bc_matrices_folder
    out:
    - compressed_folder

  compress_secondary_analysis_report_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: cellranger_multi/secondary_analysis_report_folder
    out:
    - compressed_folder

  convert_clonotypes_csv_to_tsv:
    run: ../tools/custom-bash.cwl
    in:
      input_file: cellranger_multi/clonotypes_csv
      script:
        default: |
          cat "$0" | tr "," "\t" > `basename $0 csv`tsv
    out:
    - output_file

  convert_metrics_summary_report_to_tsv:
    run: ../tools/custom-bash.cwl
    in:
      input_file: cellranger_multi/metrics_summary_report
      script:
        default: |
          cat "$0" | tr "," "\t" > `basename $0 csv`tsv
    out:
    - output_file

  cellbrowser_build:
    run: ../tools/cellbrowser-build-cellranger.cwl
    in:
      secondary_analysis_report_folder: cellranger_multi/secondary_analysis_report_folder
      filtered_feature_bc_matrix_folder: cellranger_multi/filtered_feature_bc_matrix_folder
    out:
    - html_data
    - index_html_file


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Cell Ranger Count (RNA+VDJ)"
s:name: "Cell Ranger Count (RNA+VDJ)"
s:alternateName: "Quantifies single-cell gene expression, performs V(D)J contigs assembly and clonotype calling of the sequencing data from a single 10x Genomics library in a combined manner"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/cellranger-multi.cwl
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
  Cell Ranger Count (RNA+VDJ)

  Quantifies single-cell gene expression, performs V(D)J contigs
  assembly and clonotype calling of the sequencing data from a
  single 10x Genomics library in a combined manner. The results
  of this workflow are primarily used in either “Single-Cell
  RNA-Seq Filtering Analysis”, “Single-Cell Immune Profiling Analysis”,
  or “Cell Ranger Aggregate (RNA, RNA+VDJ)” pipelines.