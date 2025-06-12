cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: DockerRequirement
  dockerPull: biowardrobe2/soupx:v0.0.1
inputs:
  feature_bc_matrices_folder:
    type: Directory
    inputBinding:
      prefix: --counts
    doc: |
      Path to the output folder produced by 'cellranger count' command
  genelist_file:
    type: File?
    inputBinding:
      prefix: --genes
    doc: |
      Path to the file with target genes (headerless, one gene per line)
  expression_threshold:
    type: float?
    inputBinding:
      prefix: --threshold
    doc: |
      Expression threshold for displaying target genes on a plot (expression > threshold)
      Default: 0.0
  fdr:
    type: float?
    inputBinding:
      prefix: --fdr
    doc: |
      FDR cutoff for expression ratio plots
      Default: 0.05
  round_counts:
    type: boolean?
    inputBinding:
      prefix: --round
    doc: |
      Round adjusted counts to integers
      Default: False
  matrix_format_version:
    type:
    - 'null'
    - type: enum
      name: matrix_format_version
      symbols:
      - '2'
      - '3'
    inputBinding:
      prefix: --format
    doc: |
      Output matrix format version. Corresponds to the latest Cell Ranger matrix format
      Default: 3
  output_prefix:
    type: string?
    inputBinding:
      prefix: --output
    doc: |
      Output prefix
      Default: soupx
outputs:
  adjusted_feature_bc_matrices_folder:
    type: Directory
    outputBinding:
      glob: '*adjusted_counts'
    doc: |
      Folder with adjusted feature-barcode matrices in MEX format
  adjusted_feature_bc_matrices_h5:
    type: File
    outputBinding:
      glob: '*adjusted_counts.h5'
    doc: |
      Adjusted feature-barcode matrices in HDF5 format
  contamination_estimation_plot:
    type: File
    outputBinding:
      glob: '*contamination_estimation_plot.pdf'
    doc: |
      Contamination estimation plot
  raw_gene_expression_plots:
    type: File?
    outputBinding:
      glob: '*raw_gene_expression_plots.pdf'
    doc: |
      Raw gene expression plots
  adjusted_gene_expression_plots:
    type: File?
    outputBinding:
      glob: '*adjusted_gene_expression_plots.pdf'
    doc: |
      Adjusted gene expression plots
  raw_gene_expression_to_pure_soup_ratio_plots:
    type: File?
    outputBinding:
      glob: '*raw_gene_expression_to_pure_soup_ratio_plots.pdf'
    doc: |
      Raw gene expression to pure soup ratio plots
  raw_to_adjusted_gene_expression_ratio_plots:
    type: File?
    outputBinding:
      glob: '*raw_to_adjusted_gene_expression_ratio_plots.pdf'
    doc: |
      Raw to adjusted gene expression ratio plots
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- run_soupx.R
stdout: soupx_stdout.log
stderr: soupx_stderr.log
label: SoupX - an R package for the estimation and removal of cell free mRNA contamination
doc: "In droplet based, single cell RNA-seq experiments, there is always a certain amount of background\nmRNAs present in the dilution that gets distributed into the droplets with cells and sequenced\nalong with them. The net effect of this is to produce a background contamination that represents\nexpression not from the cell contained within a droplet, but the solution that contained the cells.\n\nThis collection of cell free mRNAs floating in the input solution (henceforth referred to as\n“the soup”) is created from cells in the input solution being lysed. Because of this, the soup looks\ndifferent for each input solution and strongly resembles the expression pattern obtained by summing\nall the individual cells.\n\nThe aim of this package is to provide a way to estimate the composition of this soup, what fraction\nof UMIs are derived from the soup in each droplet and produce a corrected count table with the soup\nbased expression removed.\n\nThe method to do this consists of three parts:\n\n- Calculate the profile of the soup.\n- Estimate the cell specific contamination fraction.\n- Infer a corrected expression matrix.  \n"
