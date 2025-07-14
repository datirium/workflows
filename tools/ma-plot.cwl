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
      Name of column in file for the plots x-axis (ex: "baseMean")
  y_axis_column:
    type: string
    inputBinding:
      position: 7
    doc: |
      Name of column in file for the plots y-axis (ex: "log2FoldChange")
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
      glob: ./volcano_plot/MD-MA_plot
    doc: |
      Directory html data for MA-plot
  html_file:
    type: File
    outputBinding:
      glob: ./volcano_plot/MD-MA_plot/html_data/index.html
    doc: |
      HTML index file for MA-plot
  log_file_stdout:
    type: stdout
  log_file_stderr:
    type: stderr
baseCommand:
- ma_plot.sh
stdout: ma_plot-stdout.log
stderr: ma_plot-stderr.log
label: MA-plot
doc: |
  MA-plot

  Builds ma-plot from the DESeq output
