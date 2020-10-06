cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: DockerRequirement
  dockerPull: biowardrobe2/pca:v0.0.7


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
    type: string?
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
      glob: "*001.png"
    doc: "PCA1 vs PCA2 plot"

  pca2_vs_pca3_plot:
    type: File
    outputBinding:
      glob: "*002.png"
    doc: "PCA2 vs PCA3 plot"

  variance_plot:
    type: File
    outputBinding:
      glob: "*003.png"
    doc: "Variance plot"
    
  pca_3d_plot:
    type: File?
    outputBinding:
      glob: "*004.png"
    doc: "First three principal components plot"

  pca_3d_html:
    type: File?
    outputBinding:
      glob: "*.html"
    doc: "Plotly generated interactive 3D PCA plot (first three components)"

  pca_file:
    type: File
    outputBinding:
      glob: "*.tsv"
    doc: "PCA analysis results exported as TSV"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["run_pca.R"]
stderr: pca_stderr.log
stdout: pca_stdout.log