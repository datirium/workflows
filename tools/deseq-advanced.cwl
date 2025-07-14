cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap-deseq:v0.0.28
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
  fdr:
    type: float?
    inputBinding:
      prefix: --fdr
    doc: |
      In the exploratory visualization part of the analysis output only features,
      with adjusted p-value (FDR) not bigger than this value. Also the significance,
      cutoff used for optimizing the independent filtering. Default: 0.1.
  lfcthreshold:
    type: float?
    inputBinding:
      prefix: --lfcthreshold
    doc: |
      Log2 fold change threshold for determining significant differential expression.
      Genes with absolute log2 fold change greater than this threshold will be considered.
      Default: 0.59 (about 1.5 fold change)
  use_lfc_thresh:
    type: boolean
    inputBinding:
      prefix: --use_lfc_thresh
    default: true
    doc: 'Use lfcthreshold as the null hypothesis value in the results function call. Default: TRUE'
  regulation:
    type:
    - 'null'
    - type: enum
      symbols:
      - both
      - up
      - down
    inputBinding:
      prefix: --regulation
    doc: |
      Direction of differential expression comparison. β is the log2 fold change.
      - 'both' for both up and downregulated genes (|β| > lfcThreshold for greaterAbs and |β| < lfcThreshold for lessAbs, with p-values being two-tailed or maximum of the upper and lower tests, respectively).
      - 'up' for upregulated genes (β > lfcThreshold in condition2 compared to condition1).
      - 'down' for downregulated genes (β < -lfcThreshold in condition2 compared to condition1).
      Default: both
    default: both
  batchcorrection:
    type:
    - 'null'
    - type: enum
      symbols:
      - none
      - combatseq
      - limmaremovebatcheffect
    inputBinding:
      prefix: --batchcorrection
    doc: |
      Specifies the batch correction method to be applied.
      - 'combatseq' applies ComBat_seq at the beginning of the analysis, removing batch effects from the design formula before differential expression analysis.
      - 'limmaremovebatcheffect' applies removeBatchEffect from the limma package after differential expression analysis, incorporating batch effects into the model during DE analysis.
      - Default: none
    default: none
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
  threads:
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
  deseq_summary_md:
    type: File
    outputBinding:
      glob: '*summary.md'
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
baseCommand:
- run_deseq.R
stdout: deseq_stdout.log
stderr: deseq_stderr.log
doc: |
  Tool runs DESeq/DESeq2 script similar to the original one from BioWArdrobe.
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
label: deseq-advanced
