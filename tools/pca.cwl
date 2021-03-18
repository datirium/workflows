cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: DockerRequirement
  dockerPull: biowardrobe2/pca:v0.0.9


inputs:

  expression_files:
    type: File[]
    inputBinding:
      prefix: "--input"
    doc: "Input CSV/TSV files with RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand, TotalReads, Rpkm columns"

  expression_aliases:
    type:
      - "null"
      - string[]
    inputBinding:
      prefix: "--name"
    doc: "Input aliases, the order corresponds to --input order. Default: basename of --input files"

  genelist_file:
    type: File?
    inputBinding:
      prefix: "--genelist"
    doc: "Filter genes by the list from the file. Headerless, 1 gene per line"

  target_column:
    type:
    - "null"
    - type: enum
      symbols:
      - "Rpkm"
      - "TotalReads"
    inputBinding:
      prefix: "--target"
    doc: "Target column name to be used by PCA. Default: Rpkm"

  combine:
    type:
      - "null"
      - string[]
    inputBinding:
      prefix: "--combine"
    doc: "Combine inputs by columns names. Default: RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand"

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: "Output prefix. Default: pca_"


outputs:

  pca1_vs_pca2_plot:
    type: File
    outputBinding:
      glob: "*pc1_pc2_plot.png"
    doc: "PCA1 vs PCA2 plot in PNG format"
  
  pca1_vs_pca2_plot_pdf:
    type: File
    outputBinding:
      glob: "*pc1_pc2_plot.pdf"
    doc: "PCA1 vs PCA2 plot in PDF format"

  pca2_vs_pca3_plot:
    type: File
    outputBinding:
      glob: "*pc2_pc3_plot.png"
    doc: "PCA2 vs PCA3 plot in PNG format"

  pca2_vs_pca3_plot_pdf:
    type: File
    outputBinding:
      glob: "*pc2_pc3_plot.pdf"
    doc: "PCA2 vs PCA3 plot in PDF format"

  variance_plot:
    type: File
    outputBinding:
      glob: "*variance_plot.png"
    doc: "Variance plot in PNG format"
    
  variance_plot_pdf:
    type: File
    outputBinding:
      glob: "*variance_plot.pdf"
    doc: "Variance plot in PDF format"

  pca_3d_html:
    type: File?
    outputBinding:
      glob: "*pca_3d_plot.html"
    doc: "Plotly generated interactive 3D PCA plot (first three components)"

  pca_scores_file:
    type: File
    outputBinding:
      glob: "*scores.tsv"
    doc: "PCA analysis scores exported as TSV"

  pca_loadings_file:
    type: File
    outputBinding:
      glob: "*loadings.tsv"
    doc: "PCA analysis loadings exported as TSV"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["run_pca.R"]
stderr: pca_stderr.log
stdout: pca_stdout.log