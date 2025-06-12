cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: avivdemorgan/scidap-deseqspikein:v1.6.0
inputs:
  untreated_files:
    type:
    - File
    - File[]
    inputBinding:
      prefix: -u
    doc: |
      Untreated input CSV/TSV files
  treated_files:
    type:
    - File
    - File[]
    inputBinding:
      prefix: -t
    doc: |
      Treated input CSV/TSV files
  ercc_counts_treated:
    type:
    - File
    - File[]
    inputBinding:
      prefix: -ter
    doc: |
      Untreated input CSV/TSV files
  ercc_counts_untreated:
    type:
    - File
    - File[]
    inputBinding:
      prefix: -uer
    doc: |
      Untreated input CSV/TSV files
  untreated_name:
    type: string?
    inputBinding:
      prefix: -un
    doc: |
      Name for untreated condition, use only letters and numbers
  treated_name:
    type: string?
    inputBinding:
      prefix: -tn
    doc: |
      Name for treated condition, use only letters and numbers
  untreated_sample_names:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      prefix: -ua
    doc: |
      Unique aliases for untreated expression files. Default: basenames of -u without extensions
  treated_sample_names:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      prefix: -ta
    doc: |
      Unique aliases for treated expression files. Default: basenames of -t without extensions
  rpkm_cutoff:
    type: float?
    inputBinding:
      prefix: -cu
    doc: |
      Minimum threshold for rpkm filtering. Default: 5
  cluster_method:
    type:
    - 'null'
    - type: enum
      symbols:
      - row
      - column
      - both
    inputBinding:
      prefix: --cluster
    doc: |
      Hopach clustering method to be run on normalized read counts for the
      exploratory visualization part of the analysis. Default: do not run
      clustering
  row_distance:
    type:
    - 'null'
    - type: enum
      symbols:
      - cosangle
      - abscosangle
      - euclid
      - abseuclid
      - cor
      - abscor
    inputBinding:
      prefix: --rowdist
    doc: |
      Distance metric for HOPACH row clustering. Ignored if --cluster is not
      provided. Default: cosangle
  column_distance:
    type:
    - 'null'
    - type: enum
      symbols:
      - cosangle
      - abscosangle
      - euclid
      - abseuclid
      - cor
      - abscor
    inputBinding:
      prefix: --columndist
    doc: |
      Distance metric for HOPACH column clustering. Ignored if --cluster is not
      provided. Default: euclid
  center_row:
    type: boolean?
    inputBinding:
      prefix: --center
    doc: |
      Apply mean centering for feature expression prior to running
      clustering by row. Ignored when --cluster is not row or both.
      Default: do not centered
  maximum_padj:
    type: float?
    inputBinding:
      prefix: --padj
    doc: |
      In the exploratory visualization part of the analysis output only features
      with adjusted P-value not bigger than this value. Default: 0.05
  batch_file:
    type: File?
    inputBinding:
      prefix: -bf
    doc: |
      Metadata file for multi-factor analysis. Headerless TSV/CSV file.
      First column - names from --ua and --ta, second column - batch name.
      Default: None
  output_prefix:
    type: string?
    inputBinding:
      prefix: -o
    doc: |
      Output prefix. Default: deseq
  threads_count:
    type: int?
    inputBinding:
      prefix: -p
    doc: |
      Run script using multiple threads
outputs:
  diff_expr_file:
    type: File
    outputBinding:
      glob: '*report.tsv'
  read_counts_file_all:
    type: File
    outputBinding:
      glob: '*counts_all.gct'
  read_counts_file_filtered:
    type: File
    outputBinding:
      glob: '*counts_filtered.gct'
  phenotypes_file:
    type: File
    outputBinding:
      glob: '*phenotypes.cls'
  plot_lfc_vs_mean:
    type: File?
    outputBinding:
      glob: '*_ma_plot.png'
  gene_expr_heatmap:
    type: File?
    outputBinding:
      glob: '*_expression_heatmap.png'
  plot_pca:
    type: File?
    outputBinding:
      glob: '*_pca_plot.png'
  plot_lfc_vs_mean_pdf:
    type: File?
    outputBinding:
      glob: '*_ma_plot.pdf'
  gene_expr_heatmap_pdf:
    type: File?
    outputBinding:
      glob: '*_expression_heatmap.pdf'
  plot_pca_pdf:
    type: File?
    outputBinding:
      glob: '*_pca_plot.pdf'
  mds_plot_html:
    type: File?
    outputBinding:
      glob: '*_mds_plot.html'
    doc: |
      MDS plot of normalized counts. Optionally batch corrected
      based on the --remove value.
      HTML format
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
  error_msg:
    type: File?
    outputBinding:
      glob: error_msg.txt
  error_report:
    type: File?
    outputBinding:
      glob: error_report.txt
baseCommand:
- /usr/local/bin/run_deseq_for_spikein.R
stdout: deseq_stdout.log
stderr: deseq_stderr.log
doc: |
  Tool runs DESeq/DESeq2 script similar to the original one from BioWArdrobe.
  Expected input counts have already been normalized via a spike-in.
  DESeq internal normalization has been disabled in the baseCommand script call

  untreated_files and treated_files input files should have the following header (case-sensitive)
  <RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm>         - CSV
  <RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tTotalReads\tRpkm>  - TSV

  Format of the input files is identified based on file's extension
  *.csv - CSV
  *.tsv - TSV
  Otherwise used CSV by default

  The output file's rows order corresponds to the rows order of the first CSV/TSV file in
  the untreated group. Output is always saved in TSV format

  Output file includes only intersected rows from all input files. Intersected by
  RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand

  DESeq/DESeq2 always compares untreated_vs_treated groups.
  Normalized read counts and phenotype table are exported as GCT and CLS files for GSEA downstream analysis.
label: deseq-advanced-for-spikein
