cwlVersion: v1.0
class: CommandLineTool
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/visualization:v0.0.8
inputs:
  diff_expr_file:
    type: File
    inputBinding:
      position: 5
    doc: |
      TSV file holding data for the plot
  x_axis_column:
    type: string
    inputBinding:
      position: 6
    doc: |
      Name of column in file for the plots x-axis (ex: "log2FoldChange")
  y_axis_column:
    type: string
    inputBinding:
      position: 7
    doc: |
      Name of column in file for the plots y-axis (ex: "padj")
  label_column:
    type: string
    inputBinding:
      position: 8
    doc: |
      Name of column in file for each data points 'name' (ex: "GeneId")
outputs:
  html_data:
    type: Directory
    outputBinding:
      glob: ./volcano_plot/volcano_plot
    doc: |
      Directory html data for Volcano Plot
  html_file:
    type: File
    outputBinding:
      glob: ./volcano_plot/volcano_plot/html_data/index.html
    doc: |
      HTML index file for Volcano Plot
  log_file_stdout:
    type: stdout
  log_file_stderr:
    type: stderr
baseCommand:
- volcano_plot.sh
stdout: volcano_plot-stdout.log
stderr: volcano_plot-stderr.log
label: Volcano Plot
doc: |
  Volcano Plot

  Builds volcano plot from the DESeq output
